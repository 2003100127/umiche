__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import tensorflow as tf
from tensorflow.keras.utils import to_categorical
from umiche.deepconvcvae.Preprocess import preprocess as datapp
from umiche.deepconvcvae.CondiConvVAE import condiConvVAE
from umiche.plot.scatter.DimensionReduction import dimensionReduction as drplot


def sampling(
        decoder,
        num_labels,
        image_size,
        batch_size=32,
):
    xs = []
    ys = []
    # num_ = 100
    num_ = 1000
    for ilabel in range(num_labels):
        print('===>label {}'.format(ilabel))
        digit_size = image_size
        z_combo_spl = np.random.normal(loc=1, scale=4, size=[num_, 2])
        # z_combo_spl = np.random.uniform(0, 4, size=[num_, 2])
        ilabel_ = np.eye(num_labels)[np.array([ilabel])]
        for i, zcspl in enumerate(z_combo_spl):
            z_sample = np.array([zcspl])
            x_decoded = decoder.predict([z_sample, ilabel_])
            # print()
            digit = x_decoded[0].reshape(digit_size, digit_size)
            xs.append(digit.reshape(digit_size * digit_size))
            # xs.append(x_decoded[0])
            ys.append(ilabel)
            # sds = digit * 255
            # print(np.where(sds > 1, sds, 0))
    drplot().single(X=xs, y=ys, tech='TSNE', marker_size=3, cmap='tab20b', title='CondiCVAE')
    drplot().single(X=xs, y=ys, tech='PCA', marker_size=3, cmap='tab20b', title='CondiCVAE')
    drplot().single(X=xs, y=ys, tech='UMAP', marker_size=3, cmap='tab20b', title='CondiCVAE')


def data_imputation(
        image_size,
        params,
):
    pp = datapp(params['data_rv_zero_fpn'], train_ratio=3, test_ratio=2)
    (x_train, y_train), (x_test, y_test) = pp.recipeData(cv=0)
    # print(x_train.shape)
    # print(y_train.shape)
    x_train_height = x_train.shape[0]
    x_test_height = x_test.shape[0]
    x_train_width = x_train.shape[1]

    img_col_offset_size = image_size*image_size - x_train_width
    train_col_offset_mat = np.zeros([x_train_height, img_col_offset_size])
    test_col_offset_mat = np.zeros([x_test_height, img_col_offset_size])
    x_train = np.concatenate([x_train, train_col_offset_mat], axis=1)
    x_test = np.concatenate([x_test, test_col_offset_mat], axis=1)
    x_train = np.reshape(x_train, newshape=(-1, image_size, image_size))
    x_test = np.reshape(x_test, newshape=(-1, image_size, image_size))
    # print(x_train)
    # print(x_test)
    x_train = np.reshape(x_train, [-1, image_size, image_size, 1])
    x_test = np.reshape(x_test, [-1, image_size, image_size, 1])
    print(x_train.shape)
    print(y_train.shape)
    print(x_test.shape)
    print(y_test.shape)
    print(np.amax(x_train))
    print(np.amax(x_test))
    std_max_mark = np.maximum(np.amax(x_train), np.amax(x_test))
    print(std_max_mark)
    x_train = x_train.astype('float32') / std_max_mark
    x_test = x_test.astype('float32') / std_max_mark
    # print(y_train)
    # print(y_test)
    return x_train, y_train, x_test, y_test


def train(
        x_train,
        y_train,
        x_test,
        y_test,
        image_size,
        sv_fpn,
        batch_size=32,
        kernel_size=3,
        filters=16,
        latent_dim=2,
        strides=2,
        epochs=1,
):
    # num_labels = 10
    num_labels = len(np.unique(y_train))
    # print(num_labels)
    input_shape = (image_size, image_size, 1)
    label_shape = (num_labels,)

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

    outputs = decoder([encoder([inputs, y_labels])[2], y_labels])
    ccvae = tf.keras.models.Model(inputs=[inputs, y_labels], outputs=outputs, name='ccvae')
    ccvae.add_loss(condicvae.loss(outputs, z_mean, z_log_var))
    ccvae.compile(optimizer='rmsprop')
    ccvae.summary()

    ccvae.fit(
        [x_train, to_categorical(y_train)],
        epochs=epochs,
        batch_size=batch_size,
        validation_data=([x_test, to_categorical(y_test)], None),
    )

    # model_name = "cvae_model20"
    # save_dir = "cvae_model20"
    # if not os.path.isdir(save_dir):
    #     os.makedirs(save_dir)
    # filename = model_name + '.tf'
    # filepath = os.path.join(save_dir, filename)
    # ccvae.save_weights(filepath)

    ccvae.save(sv_fpn)

    sampling(
        decoder,
        num_labels=num_labels,
        batch_size=batch_size,
        image_size=image_size,
    )
    return


if __name__ == "__main__":
    from umiche.path import to

    params = {
        'data_fpn': to('data/deepconvcvae/flt.h5'),
        'data_rv_zero_fpn': to('data/deepconvcvae/flt_rv_zero_col.h5'),
    }

    image_size = 144

    x_train, y_train, x_test, y_test = data_imputation(
        params=params,
        image_size=image_size,
    )

    print(train(
        x_train=x_train,
        y_train=y_train,
        x_test=x_test,
        y_test=y_test,
        image_size=image_size,
        batch_size=32,
        kernel_size=3,
        filters=16,
        latent_dim=2,
        strides=2,
        epochs=1,
        sv_fpn='./cvae',
    ))