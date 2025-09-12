__version__ = "v1.3"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Dict, Optional, Tuple, List

import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse import issparse
import tensorflow as tf
from tensorflow.keras.utils import to_categorical
from tensorflow.keras import backend as K

from umiche.deepconvcvae.CondiConvVAE import condiConvVAE


def _to_dense(X):
    return X.A if hasattr(X, "A") else (X.toarray() if issparse(X) else X)


class DeepConvCVAEAdapter:
    def __init__(
        self,
        image_size: Optional[int] = None,
        gene_strategy: str = "hvg",
        max_genes_cap: Optional[int] = 20000,
        min_gene_counts: int = 1,
        align_multiple: int = 4,
        latent_dim: int = 2,
        filters: int = 16,
        kernel_size: int = 3,
        strides: int = 2,
        epochs: int = 50,
        batch_size: int = 64,
        optimizer: str = "rmsprop",
        z_loc: float = 0.0,
        z_scale: float = 1.0,
        sim_mode: str = "zinb",
        use_posterior_z: bool = True,
        auto_tau: bool = True,
        tau: Optional[float] = None,
        tau_grid: Tuple[float, ...] = (0.6, 0.8, 1.0, 1.2, 1.6, 2.0),
        tau_pilot_n: int = 200,
        theta_min: float = 1.0,
        theta_max: float = 100.0,
        zi_cap: float = 0.98,
        gamma: float = 0.9,
        alpha_scale: float = 0.035,
        theta_cell_sigma: float = 0.35,
        sigma_cell: float = 0.25,
        sigma_gene: float = 0.10,
        seed: int = 1,
        verbose: int = 1,
    ):
        self.image_size = image_size
        self.gene_strategy = gene_strategy
        self.max_genes_cap = max_genes_cap
        self.min_gene_counts = min_gene_counts
        self.align_multiple = max(1, int(align_multiple))

        self.latent_dim = int(latent_dim)
        self.filters = int(filters)
        self.kernel_size = int(kernel_size)
        self.strides = int(strides)
        self.epochs = int(epochs)
        self.batch_size = int(batch_size)
        self.optimizer = optimizer

        self.z_loc = float(z_loc)
        self.z_scale = float(z_scale)

        self.sim_mode = sim_mode
        self.use_posterior_z = bool(use_posterior_z)

        self.auto_tau = bool(auto_tau)
        self.tau = tau
        self.tau_grid = tuple(tau_grid)
        self.tau_pilot_n = int(tau_pilot_n)

        self.theta_min = float(theta_min)
        self.theta_max = float(theta_max)
        self.zi_cap = float(zi_cap)

        self.gamma = float(gamma)
        self.alpha_scale = float(alpha_scale)
        self.theta_cell_sigma = float(theta_cell_sigma)
        self.sigma_cell = float(sigma_cell)
        self.sigma_gene = float(sigma_gene)

        self.seed = int(seed)
        self.verbose = int(verbose)

        self.square_genes_: Optional[np.ndarray] = None
        self.genes_: Optional[np.ndarray] = None
        self.scale_max_: Optional[float] = None
        self.label_levels_: Optional[np.ndarray] = None
        self.label_to_index_: Optional[Dict[str, int]] = None
        self.num_labels_: Optional[int] = None

        self.encoder_ = None
        self.decoder_ = None
        self.decoder_gen_ = None
        self.model_ = None

        self._libs_by_label: Dict[str, np.ndarray] = {}
        self._czero_median_train: float = 0.0
        self._z_bank_mu: Dict[str, np.ndarray] = {}
        self._z_bank_std: Dict[str, np.ndarray] = {}

        self.nb_mu_: Dict[str, np.ndarray] = {}
        self.nb_theta_: Dict[str, np.ndarray] = {}
        self.nb_pi_: Dict[str, np.ndarray] = {}

    def _choose_square_genes(self, adata: ad.AnnData, strategy: str,
                             cap: Optional[int], min_counts: int) -> Tuple[np.ndarray, int]:
        at = adata.copy()
        sc.pp.filter_genes(at, min_counts=min_counts)
        if at.n_vars == 0:
            raise ValueError("No genes left after filtering.")
        M = min(cap, at.n_vars) if cap is not None else at.n_vars
        s = int(np.floor(np.sqrt(M)))
        s = max(self.align_multiple, (s // self.align_multiple) * self.align_multiple)
        n_select = s * s
        if strategy == "hvg":
            sc.pp.highly_variable_genes(at, n_top_genes=min(n_select, at.n_vars),
                                        flavor="seurat_v3", subset=True)
            if at.n_vars < n_select:
                X = _to_dense(at.X); m = X.mean(0)
                idx = np.argsort(m)[::-1][:n_select]
                genes = at.var_names.values[idx]
            else:
                genes = at.var_names.values
        else:
            X = _to_dense(at.X); m = X.mean(0)
            idx = np.argsort(m)[::-1][:n_select]
            genes = at.var_names.values[idx]
        return genes, s

    def _prepare_images_and_labels(self, adata: ad.AnnData, label_col: str,
                                   genes: np.ndarray, for_fit: bool):
        X = _to_dense(adata[:, genes].X).astype(np.float32, copy=False)
        X[X < 0] = 0
        if for_fit:
            self.scale_max_ = float(np.max(X)) if np.max(X) > 0 else 1.0
        scale_max = self.scale_max_ if self.scale_max_ else 1.0
        X = X / scale_max
        H = W = self.image_size
        X = X.reshape(-1, H, W, 1)

        y_raw = adata.obs[label_col].astype(str).values
        if for_fit:
            self.label_levels_ = np.unique(y_raw)
            self.label_to_index_ = {lab: i for i, lab in enumerate(self.label_levels_)}
            self.num_labels_ = len(self.label_levels_)
        idx = np.array([self.label_to_index_[lab] for lab in y_raw], dtype=np.int32)
        return X, idx

    def fit(self, adata_train: ad.AnnData, label_col: str = "cell_type",
            genes: Optional[np.ndarray] = None, val_data=None):
        tf.keras.utils.set_random_seed(self.seed)

        if genes is not None:
            genes = np.array(genes, dtype=str)
            s = int(np.floor(np.sqrt(len(genes))))
            s = max(self.align_multiple, (s // self.align_multiple) * self.align_multiple)
            genes = genes[: s * s]
            self.image_size = s
            self.square_genes_ = genes
        else:
            genes, s = self._choose_square_genes(adata_train, self.gene_strategy,
                                                 self.max_genes_cap, self.min_gene_counts)
            self.square_genes_ = genes
            self.image_size = s
        self.genes_ = self.square_genes_

        X_tr, y_tr = self._prepare_images_and_labels(adata_train, label_col, self.genes_, for_fit=True)
        y1h_tr = to_categorical(y_tr, num_classes=len(self.label_levels_))

        input_shape = (self.image_size, self.image_size, 1)
        label_shape = (len(self.label_levels_),)
        inputs = tf.keras.layers.Input(shape=input_shape, name="encoder_input")
        y_labels = tf.keras.layers.Input(shape=label_shape, name="class_labels")

        cvae = condiConvVAE(
            image_size=self.image_size,
            input_shape=input_shape,
            label_shape=label_shape,
            batch_size=self.batch_size,
            kernel_size=self.kernel_size,
            filters=self.filters,
            latent_dim=self.latent_dim,
            strides=self.strides,
            epochs=self.epochs,
            inputs=inputs,
            y_labels=y_labels,
        )
        shape, z_mean, z_log_var, encoder = cvae.encoding()
        decoder = cvae.decoding(shape)

        z_post = encoder([inputs, y_labels])[2]
        outputs = decoder([z_post, y_labels])
        in_h = in_w = self.image_size
        oh, ow = K.int_shape(outputs)[1:3]
        if (oh is not None and oh != in_h) or (ow is not None and ow != in_w):
            dh = (oh - in_h) if oh is not None else 0
            dw = (ow - in_w) if ow is not None else 0
            if (dh > 0) or (dw > 0):
                outputs = tf.keras.layers.Cropping2D(((0, max(dh, 0)), (0, max(dw, 0))))(outputs)
            if (dh < 0) or (dw < 0):
                outputs = tf.keras.layers.ZeroPadding2D(((0, max(-dh, 0)), (0, max(-dw, 0))))(outputs)

        model = tf.keras.models.Model([inputs, y_labels], outputs, name="ccvae")
        model.add_loss(cvae.loss(outputs, z_mean, z_log_var))
        model.compile(optimizer=self.optimizer)
        if self.verbose:
            model.summary()
        model.fit([X_tr, y1h_tr], epochs=self.epochs, batch_size=self.batch_size,
                  validation_data=None if val_data is None else val_data, verbose=self.verbose)

        self.encoder_ = encoder
        self.decoder_ = decoder
        self.model_ = model

        z_in = tf.keras.layers.Input(shape=(self.latent_dim,), name="z_infer")
        y_in = tf.keras.layers.Input(shape=(len(self.label_levels_),), name="y_infer")
        dec_o = decoder([z_in, y_in])
        oh2, ow2 = K.int_shape(dec_o)[1:3]
        if (oh2 is not None and oh2 != in_h) or (ow2 is not None and ow2 != in_w):
            dh = (oh2 - in_h) if oh2 is not None else 0
            dw = (ow2 - in_w) if ow2 is not None else 0
            if (dh > 0) or (dw > 0):
                dec_o = tf.keras.layers.Cropping2D(((0, max(dh, 0)), (0, max(dw, 0))))(dec_o)
            if (dh < 0) or (dw < 0):
                dec_o = tf.keras.layers.ZeroPadding2D(((0, max(-dh, 0)), (0, max(-dw, 0))))(dec_o)
        self.decoder_gen_ = tf.keras.models.Model([z_in, y_in], dec_o, name="decoder_gen")

        X_counts = _to_dense(adata_train[:, self.genes_].X)
        X_counts = np.asarray(X_counts)
        libs = X_counts.sum(axis=1).astype(np.int64).ravel()
        labels = adata_train.obs[label_col].astype(str).values
        for lab in self.label_levels_:
            sel = (labels == lab)
            vals = libs[sel]
            if vals.size == 0:
                vals = libs
            self._libs_by_label[lab] = vals
        cz_vec = np.mean(X_counts == 0, axis=1)
        self._czero_median_train = float(np.median(cz_vec))

        if self.use_posterior_z:
            mu, logvar, _ = encoder([X_tr, y1h_tr])
            mu = mu.numpy(); std = np.sqrt(np.exp(logvar.numpy()))
            for i, lab in enumerate(self.label_levels_):
                idx = np.where(y_tr == i)[0]
                self._z_bank_mu[lab] = mu[idx]
                self._z_bank_std[lab] = std[idx]

        for lab in self.label_levels_:
            C = X_counts[labels == lab]
            if C.shape[0] == 0:
                C = X_counts
            L = C.sum(1)
            L_med = np.median(L) if np.median(L) > 0 else 1.0
            sf = np.maximum(L / L_med, 1e-8)
            Y = C / sf[:, None]
            mu_g = Y.mean(0) + 1e-8
            var_g = Y.var(0) + 1e-8
            theta_g = (mu_g**2) / np.maximum(var_g - mu_g, 1e-8)
            theta_g = np.clip(theta_g, self.theta_min, self.theta_max)
            z_emp = (C == 0).mean(0)
            z_nb = (theta_g / (theta_g + mu_g)) ** theta_g
            pi_g = (z_emp - z_nb) / np.maximum(1.0 - z_nb, 1e-8)
            pi_g = np.clip(pi_g, 0.0, self.zi_cap)

            self.nb_mu_[lab] = mu_g.astype(np.float32)
            self.nb_theta_[lab] = theta_g.astype(np.float32)
            self.nb_pi_[lab] = pi_g.astype(np.float32)

        return self

    def _sample_z(self, lab: str, k: int) -> np.ndarray:
        if self.use_posterior_z and lab in self._z_bank_mu and self._z_bank_mu[lab].size:
            idx = np.random.randint(0, self._z_bank_mu[lab].shape[0], size=k)
            mu = self._z_bank_mu[lab][idx]
            std = self._z_bank_std[lab][idx]
            eps = np.random.normal(size=(k, self.latent_dim)).astype(np.float32)
            return (mu + std * eps).astype(np.float32)
        return np.random.normal(loc=self.z_loc, scale=self.z_scale,
                                size=(k, self.latent_dim)).astype(np.float32)

    @staticmethod
    def _softmax_temp(X: np.ndarray, tau: float, axis: int = -1, eps: float = 1e-8) -> np.ndarray:
        X = X / max(1e-6, tau)
        X = X - X.max(axis=axis, keepdims=True)
        np.exp(X, out=X)
        s = X.sum(axis=axis, keepdims=True) + eps
        return X / s

    def _auto_tune_tau(self) -> float:
        n = self.tau_pilot_n
        per_lab = max(1, n // len(self.label_levels_))
        H = W = self.image_size
        G = H * W
        cz_target = self._czero_median_train
        all_scores = []
        for tau in self.tau_grid:
            cz = []
            for i, lab in enumerate(self.label_levels_):
                k = per_lab
                z = self._sample_z(lab, k)
                y1h = np.zeros((k, len(self.label_levels_)), dtype=np.float32); y1h[:, i] = 1.0
                V = self.decoder_gen_.predict([z, y1h], batch_size=256, verbose=0).reshape(k, G)
                P = self._softmax_temp(V, tau=tau, axis=1)
                L = int(np.median(self._libs_by_label[lab]))
                Xs = np.vstack([np.random.multinomial(L, P[j]) for j in range(k)])
                cz.append((Xs == 0).mean())
            all_scores.append((abs(np.median(cz) - cz_target), tau))
        best_tau = min(all_scores, key=lambda x: x[0])[1]
        if self.verbose:
            print(f"[DeepConvCVAE] auto tau -> {best_tau:.3f}")
        return best_tau

    def simulate(self,
                 target_counts: Optional[Dict[str, int]] = None,
                 n_cell_new: Optional[int] = None,
                 label_col_out: str = "cell_type",
                 batch_size: int = 256) -> ad.AnnData:
        assert self.decoder_gen_ is not None and self.genes_ is not None and self.num_labels_ is not None, \
            "Please call fit() first."
        tf.keras.utils.set_random_seed(self.seed)

        if target_counts is not None:
            order = list(self.label_levels_)
            per_label = [int(target_counts.get(lab, 0)) for lab in order]
        elif n_cell_new is not None:
            base = n_cell_new // len(self.label_levels_)
            per_label = [base] * len(self.label_levels_)
            per_label[0] += (n_cell_new - base * len(self.label_levels_))
        else:
            raise ValueError("Provide target_counts or n_cell_new.")

        H = W = self.image_size
        G = H * W
        X_list: List[np.ndarray] = []
        y_list: List[str] = []

        if self.sim_mode == "multinomial" and (self.tau is None) and self.auto_tau:
            self.tau = self._auto_tune_tau()
        tau = 1.0 if self.tau is None else float(self.tau)

        for i, lab in enumerate(self.label_levels_):
            k = per_label[i]
            if k <= 0:
                continue

            z = self._sample_z(lab, k)
            y1h = np.zeros((k, len(self.label_levels_)), dtype=np.float32); y1h[:, i] = 1.0
            V = self.decoder_gen_.predict([z, y1h], batch_size=batch_size, verbose=0).reshape(k, G)

            libs = self._libs_by_label.get(lab, None)
            if libs is None or libs.size == 0:
                libs = np.array([G], dtype=int)
            Ls = np.random.choice(libs, size=k, replace=True).astype(np.int64)
            Ls[Ls < 1] = 1

            if self.sim_mode == "multinomial":
                P = self._softmax_temp(V, tau=tau, axis=1)
                Xs = np.vstack([np.random.multinomial(int(Ls[j]), P[j]) for j in range(k)]).astype(np.int64)

            else:
                mu_g = self.nb_mu_[lab]  # (G,)
                th_g = self.nb_theta_[lab]  # (G,)
                pi_g = self.nb_pi_[lab]  # (G,)

                w = np.maximum(V, 0)
                w_sum = w.sum(axis=1, keepdims=True)
                bad = (w_sum.squeeze() == 0)
                if np.any(bad):
                    w[bad, :] = mu_g[None, :]
                    w_sum[bad, :] = w[bad, :].sum(axis=1, keepdims=True)
                Pw = w / w_sum
                mu_target = (Ls.reshape(-1, 1) * Pw).astype(np.float32)

                L_med = np.median(libs) if np.median(libs) > 0 else 1.0
                sf = (Ls / L_med).reshape(-1, 1).astype(np.float32)
                mu_base = sf * mu_g.reshape(1, -1)
                mu = 0.7 * mu_target + 0.3 * mu_base
                mu = np.maximum(mu, 1e-8)

                if self.sigma_cell > 0:
                    cell_noise = np.exp(np.random.normal(0.0, self.sigma_cell, size=(k, 1))).astype(np.float32)
                    mu *= cell_noise
                if self.sigma_gene > 0:
                    gene_noise = np.exp(np.random.normal(0.0, self.sigma_gene, size=(1, G))).astype(np.float32)
                    mu *= gene_noise
                mu *= (Ls.reshape(-1, 1) / np.maximum(mu.sum(axis=1, keepdims=True), 1e-8))

                if self.theta_cell_sigma > 0:
                    kdisp = np.exp(np.random.normal(0.0, self.theta_cell_sigma, size=(k, 1))).astype(np.float32)
                else:
                    kdisp = np.ones((k, 1), dtype=np.float32)

                out_rows = []
                step = max(1, 8192 // G)
                for st in range(0, k, step):
                    ed = min(k, st + step)
                    mu_sub = mu[st:ed]
                    th_eff = th_g.reshape(1, -1) / kdisp[st:ed]
                    th_eff = np.clip(th_eff, 0.05, self.theta_max)

                    rate = np.random.gamma(shape=th_eff, scale=mu_sub / th_eff)
                    Xp = np.random.poisson(rate).astype(np.int64)

                    L_ratio = (L_med / Ls[st:ed]).reshape(-1, 1)
                    pi_eff = 1.0 - np.power(1.0 - pi_g.reshape(1, -1), np.power(L_ratio, self.gamma))
                    pi_eff = np.clip(pi_eff, 0.0, 0.995)
                    drop_mask = (np.random.rand(ed - st, G) < pi_eff)
                    Xp[drop_mask] = 0

                    need = (Ls[st:ed].reshape(-1, 1) - Xp.sum(axis=1, keepdims=True)).astype(np.int64)
                    need = np.maximum(need, 0).ravel()
                    if np.any(need > 0):
                        base = mu_sub + 1e-8
                        base /= base.sum(axis=1, keepdims=True)
                        for jj, add_n in enumerate(need):
                            if add_n <= 0: continue
                            alpha_total = max(1.0, self.alpha_scale * float(Ls[st + jj]))
                            alpha = base[jj] * alpha_total
                            dirv = np.random.gamma(shape=alpha, scale=1.0)
                            p = dirv / np.sum(dirv)
                            add = np.random.multinomial(int(add_n), p)
                            Xp[jj] += add

                    out_rows.append(Xp)

                Xs = np.vstack(out_rows)

            X_list.append(Xs)
            y_list.extend([lab] * k)

        X_counts = np.vstack(X_list) if X_list else np.zeros((0, G), dtype=np.int64)
        adata_sim = ad.AnnData(X_counts)
        adata_sim.var_names = self.genes_.astype(str)
        adata_sim.obs[label_col_out] = np.array(y_list, dtype=object)
        return adata_sim

    def fit_and_simulate(self, adata_train: ad.AnnData, label_col: str = "cell_type",
                         target_counts: Optional[Dict[str, int]] = None,
                         n_cell_new: Optional[int] = None, val_data=None) -> ad.AnnData:
        self.fit(adata_train, label_col=label_col, val_data=val_data)
        return self.simulate(target_counts=target_counts, n_cell_new=n_cell_new, label_col_out=label_col)

    @staticmethod
    def counts_from_reference(adata_ref: ad.AnnData, label_col: str = "cell_type") -> Dict[str, int]:
        vc = adata_ref.obs[label_col].astype(str).value_counts()
        return {k: int(v) for k, v in vc.items()}
