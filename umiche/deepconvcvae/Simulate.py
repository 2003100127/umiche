__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import tensorflow as tf
from umiche.deepconvcvae.CondiConvVAE import condiConvVAE
from umiche.plot.scatter.DimensionReduction import dimensionReduction as drplot


def sampling(
        decoder,
        num_labels,
        model_name,
        image_size,
        num_dots=1000,
        batch_size=32,
        sv_fpn=None,
):
    xs = []
    ys = []
    for ilabel in range(num_labels):
        print('===>label {}'.format(ilabel))
        digit_size = image_size
        # z_combo_spl = np.random.poisson(5, size=[num_dots, 2])
        z_combo_spl = np.random.normal(loc=1, scale=2, size=[num_dots, 2])
        # z_combo_spl = np.random.normal(loc=1, scale=4, size=[num_dots, 2])
        # z_combo_spl = np.random.uniform(0, 4, size=[num_dots, 2])
        ilabel_ = np.eye(num_labels)[np.array([ilabel])]
        for i, zcspl in enumerate(z_combo_spl):
            z_sample = np.array([zcspl])
            x_decoded = decoder.predict([z_sample, ilabel_])
            # x_decoded = decoder.get_layer("decoder").predict([z_sample, ilabel_])
            # print(x_decoded.shape)
            # print(x_decoded[0].shape)
            digit = x_decoded[0].reshape(digit_size, digit_size)
            # print(digit)
            t = digit.reshape(digit_size * digit_size)*260
            t = np.floor(t).astype(int)
            # print(t[t>1])
            # print(t[t>1].shape)

            xs.append(t)
            # xs.append(x_decoded[0])
            ys.append(ilabel)
            # sds = digit * 255
            # print(np.where(sds > 1, sds, 0))
    print(np.array(xs))
    from scipy.sparse import csr_matrix

    print(csr_matrix(np.array(xs)))
    import pandas as pd
    sdf = pd.DataFrame(np.array(xs))
    sdf['Y'] = ys
    # sdf = pd.DataFrame.sparse.from_spmatrix(csr_matrix(np.array(xs)))
    print(sdf)
    if sv_fpn:
        # sv_h5_fpn = 'data_' + str(num_dots) + '_' + model_name + '.h5'
        sdf.to_hdf(sv_fpn, key='df', mode='w')
        print('saved')
        # print(pd.read_hdf(sv_h5_fpn, 'df'))

    # import pandas as pd
    # df = pd.read_hdf('data.h5', 'df')
    # df = df.drop(columns=['Y'])
    # import numpy as np
    # print(np.unique(df.to_numpy().flatten()))

    # drplot().single(X=xs, y=ys, tech='TSNE', marker_size=3, cmap='tab20b', title='CondiCVAE')
    # drplot().single(X=xs, y=ys, tech='PCA', marker_size=3, cmap='tab20b', title='CondiCVAE')
    # drplot().single(X=xs, y=ys, tech='UMAP', marker_size=3, cmap='tab20b', title='CondiCVAE')


def simulate(
        model_fpn,
        image_size,
        batch_size=32,
        kernel_size=3,
        filters=16,
        latent_dim=2,
        strides=2,
        epochs=1,
        num_dots=100,
        model_name='cvae_model50',
        sv_fpn=None,
):
    # model_name = "./cvae_model20"
    # save_dir = "./cvae_model20"
    # if not os.path.isdir(save_dir):
    #     os.makedirs(save_dir)
    # filename = model_name + '.tf'
    # filepath = os.path.join(save_dir, filename)
    # # ccvae.save_weights(filepath)
    # print(filepath)
    num_labels = 11
    input_shape = (image_size, image_size, 1)
    label_shape = (num_labels, )

    inputs = tf.keras.layers.Input(shape=input_shape, name='encoder_input')
    y_labels = tf.keras.layers.Input(shape=label_shape, name='class_labels')
    print(inputs)
    print(y_labels)
    condicvae = condiConvVAE(
        image_size=image_size,
        input_shape=input_shape,
        label_shape=label_shape,
        batch_size=batch_size,
        kernel_size=kernel_size,
        filters=filters,
        latent_dim=latent_dim,
        strides=strides,
        epochs=epochs,
        inputs=inputs,
        y_labels=y_labels,
    )
    shape, z_mean, z_log_var, encoder = condicvae.encoding()
    # encoder.summary()
    decoder = condicvae.decoding(shape)
    # decoder.summary()
    # plot_model(decoder, to_file='cvae_cnn_decoder.png', show_shapes=True)

    outputs = decoder([encoder([inputs, y_labels])[2], y_labels])
    ccvae = tf.keras.models.Model(inputs=[inputs, y_labels], outputs=outputs, name='ccvae')

    ccvae.load_weights(
        model_fpn,
    )

    sampling(
        decoder,
        num_labels=num_labels,
        model_name=model_name,
        image_size=image_size,
        num_dots=num_dots,
        batch_size=32,
        sv_fpn=sv_fpn,
    )


if __name__ == "__main__":
    print(simulate(
        image_size=144,
        batch_size=32,
        kernel_size=3,
        filters=16,
        latent_dim=2,
        num_dots=100,
        strides=2,
        epochs=1,
        model_fpn='../data/deepconvcvae/cvae_model50/cvae_model50.tf',
        sv_fpn='../data/deepconvcvae/cvae_model50/cvae_model50.h5',
    ))