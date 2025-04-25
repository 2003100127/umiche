__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from umiche.deepconvcvae.Simulate import simulate
from umiche.deepconvcvae.Train import train, data_imputation


def data_in(
        params,
        image_size,
):
    return data_imputation(
        image_size=image_size,
        params=params,
    )


def train_deepconvcvae(
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
    return train(
        x_train=x_train,
        y_train=y_train,
        x_test=x_test,
        y_test=y_test,
        image_size=image_size,
        batch_size=batch_size,
        kernel_size=kernel_size,
        filters=filters,
        latent_dim=latent_dim,
        strides=strides,
        epochs=epochs,
        sv_fpn=sv_fpn,
    )


def simu_deepconvcvae(
        image_size,
        batch_size,
        kernel_size,
        filters,
        latent_dim,
        strides,
        epochs,
        num_dots,
        model_fpn,
        sv_fpn,
):
    return simulate(
        image_size=image_size,
        batch_size=batch_size,
        kernel_size=kernel_size,
        filters=filters,
        num_dots=num_dots,
        latent_dim=latent_dim,
        strides=strides,
        epochs=epochs,
        model_fpn=model_fpn,
        sv_fpn=sv_fpn,
    )


if __name__ == '__main__':
    from umiche.path import to

    train_deepconvcvae()

    # print(simu_deepconvcvae(
    #     image_size=144,
    #     batch_size=32,
    #     kernel_size=3,
    #     filters=16,
    #     latent_dim=2,
    #     strides=2,
    #     epochs=1,
    #     num_dots=100,
    #     model_fpn='../data/deepconvcvae/cvae_model50/cvae_model50.tf',
    #     sv_fpn='../data/deepconvcvae/cvae_model50/cvae_model50.h5',
    # ))

    x_train, y_train, x_test, y_test = data_in(
        params={
            'data_fpn': to('data/deepconvcvae/flt.h5'),
            'data_rv_zero_fpn': to('data/deepconvcvae/flt_rv_zero_col.h5'),
        },
        image_size=144,
    )

    print(train_deepconvcvae(
        x_train=x_train,
        y_train=y_train,
        x_test=x_test,
        y_test=y_test,
        image_size=144,
        batch_size=32,
        kernel_size=3,
        filters=16,
        latent_dim=2,
        strides=2,
        epochs=1,
        sv_fpn='./cvae',
    ))