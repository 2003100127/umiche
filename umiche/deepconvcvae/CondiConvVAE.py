__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"


import tensorflow as tf
from tensorflow.keras import backend as K


class condiConvVAE:

    def __init__(
            self,
            input_shape,
            image_size,
            label_shape,
            batch_size,
            kernel_size,
            filters,
            latent_dim,
            strides,
            epochs,
            inputs,
            y_labels,
    ):
        self.image_size = image_size
        self.input_shape = input_shape
        self.label_shape = label_shape
        self.batch_size = batch_size
        self.kernel_size = kernel_size
        self.filters = filters
        self.latent_dim = latent_dim
        self.strides = strides
        self.epochs = epochs
        self.inputs = inputs
        self.y_labels = y_labels

    def reparameterize(self, params):
        z_mean, z_log_var = params
        batch = K.shape(z_mean)[0]
        dim = K.int_shape(z_mean)[1]
        epsilon = K.random_normal(shape=(batch, dim))
        return z_mean + K.exp(0.5 * z_log_var) * epsilon

    def loss(self, outputs, z_mean, z_log_var):
        # VAE loss = mse_loss or xent_loss + kl_loss
        beta = 1.0
        # rloss = tf.keras.losses.binary_crossentropy(K.flatten(inputs), K.flatten(outputs))
        rloss = tf.keras.losses.mse(K.flatten(self.inputs), K.flatten(outputs))
        rloss *= self.image_size * self.image_size
        kl_loss = 1 + z_log_var - K.square(z_mean) - K.exp(z_log_var)
        kl_loss = K.sum(kl_loss, axis=-1)
        kl_loss *= -0.5 * beta
        return K.mean(rloss + kl_loss)

    def encoding(self, ):
        x = tf.keras.layers.Dense(self.image_size * self.image_size)(self.y_labels)
        x = tf.keras.layers.Reshape((self.image_size, self.image_size, 1))(x)
        x = tf.keras.layers.concatenate([self.inputs, x])
        x = tf.keras.layers.Conv2D(self.filters, self.kernel_size, self.strides, padding='same', activation='relu')(x)
        x = tf.keras.layers.Conv2D(self.filters * 2, self.kernel_size, self.strides, padding='same', activation='relu')(x)
        # for decoder
        shape = K.int_shape(x)
        # generate latent vector Q(z|X)
        x = tf.keras.layers.Flatten()(x)
        x = tf.keras.layers.Dense(16, activation='relu')(x)
        z_mean = tf.keras.layers.Dense(self.latent_dim, name='z_mean')(x)
        z_log_var = tf.keras.layers.Dense(self.latent_dim, name='z_log_var')(x)

        # use reparameterization trick to push the sampling out as input
        z = tf.keras.layers.Lambda(self.reparameterize, output_shape=(self.latent_dim,), name='z',)([z_mean, z_log_var])
        encoder = tf.keras.models.Model(inputs=[self.inputs, self.y_labels], outputs=[z_mean, z_log_var, z], name='encoder')
        return shape, z_mean, z_log_var, encoder

    def decoding(self, shape,):
        latent_inputs = tf.keras.layers.Input(shape=(self.latent_dim,), name='z_sampling')
        x = tf.keras.layers.concatenate([latent_inputs, self.y_labels])
        x = tf.keras.layers.Dense(shape[1]*shape[2]*shape[3], activation='relu')(x)
        x = tf.keras.layers.Reshape((shape[1], shape[2], shape[3]))(x)
        x = tf.keras.layers.Conv2DTranspose(self.filters, self.kernel_size, self.strides, padding='same', activation='relu')(x)
        x = tf.keras.layers.Conv2DTranspose(self.filters / 2, self.kernel_size, self.strides, padding='same', activation='relu')(x)
        x = tf.keras.layers.Conv2DTranspose(1, self.kernel_size, padding='same', activation='sigmoid', name='decoder_output')(x)
        decoder = tf.keras.models.Model(inputs=[latent_inputs, self.y_labels], outputs=x, name='decoder')
        return decoder