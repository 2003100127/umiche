# VAE

A variational autoencoder (VAE)[^1] is a generative model in machine learning designed to encode input data into a latent space and subsequently decode it to reconstruct the original input. VAEs excel at generating new data samples that resemble a given dataset and are also valuable for various unsupervised learning tasks.

[^1]: Cemgil, T., Ghaisas, S., Dvijotham, K., Gowal, S., & Kohli, P. (2020). The autoencoding variational autoencoder. Advances in Neural Information Processing Systems, 33, 15077-15087.

We utilised a variant of the variational autoencoder (VAE), the conditional VAE (CVAE)[^2], to simulate transcriptomics data conditioned on cell types. The CVAE allows for training simulators on highly dimensional, yet sparse single-cell gene expression data at a significantly accelerated speed. Additionally, it integrates various types of neural networks into the parameterized statistical inference process, offering greater flexibility than purely statistical inference methods for denoising multimodally distributed data.

[^2]: Lopez-Martin, M., Carro, B., Sanchez-Esguevillas, A., & Lloret, J. (2017). Conditional variational autoencoder for prediction and feature recovery applied to intrusion detection in iot. Sensors, 17(9), 1967.

!!! abstract "Background: scRNA-seq count matrix simulation"

    The recent surge of interest in scRNA-seq stems from its ability to elucidate the correlation between genetics and diseases at the level of individual cell types with greater granularity. To discern the cell-specific gene expression landscapes, numerous scRNA-seq analysis tools have emerged. The reliability of these tools in practical applications hinges on their evaluation against realistically labelled data. Consequently, many simulators have been developed to ensure the availability of ground truth data. Most of these simulations are based on statistical inference techniques, using well-estimated parameters of probability distributions to characterise given scRNA-seq data. However, statistical inference methods have two limitations: long training time and a lack of denoising strategies. Firstly, methods in this category typically lack technical support from computing acceleration technologies, such as graphics processing units (GPUs), particularly when compared to deep learning computing libraries. However, the practical utility of simulated data varies across cell, tissue, and disease types, potentially necessitating model retraining to accommodate a range of application conditions. On the other hand, the inherent variability of complex biological systems often eludes capture by rigid statistical inference methods that rely solely on fixed probability distribution assumptions. This limitation may constrain their ability to effectively denoise data exhibiting multimodal distributions.

!!! danger "Point"

    A count matrix can be simulated on-demand.

# Programming CVAE

We implemented a CVAE framework based on convolutional neural networks.

!!! warning

    Please note that before using it, we need to install [`tensorflow`](https://www.tensorflow.org).:material-arrow-down-right:

    ``` shell
    pip install tensorflow
    ```

The following vignette shows the framework.

:material-language-python: `Python`
``` py linenums="1"
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
```

# Simulation

We have trained a model which we can use to simulate 100 cells (`num_dots=100`) per cell cluster. The final simulated count matrix will be saved in `h5` format.

:material-language-python: `Python`
``` py linenums="1"
import umiche as uc

uc.mat.simu_deepconvcvae(
    image_size=144,
    batch_size=32,
    kernel_size=3,
    filters=16,
    latent_dim=2,
    strides=2,
    epochs=1,
    num_dots=100,
    model_fpn='../data/deepconvcvae/cvae_model50/cvae_model50.tf',
    sv_fpn='../data/deepconvcvae/cvae_model50/cvae_model50.h5',
)
```

:material-console: `console`
``` shell
KerasTensor(type_spec=TensorSpec(shape=(None, 144, 144, 1), dtype=tf.float32, name='encoder_input'), name='encoder_input', description="created by layer 'encoder_input'")
KerasTensor(type_spec=TensorSpec(shape=(None, 11), dtype=tf.float32, name='class_labels'), name='class_labels', description="created by layer 'class_labels'")
2024-07-30 22:41:02.531703: I tensorflow/core/platform/cpu_feature_guard.cc:182] This TensorFlow binary is optimized to use available CPU instructions in performance-critical operations.
To enable the following instructions: SSE SSE2 SSE3 SSE4.1 SSE4.2 AVX AVX2 FMA, in other operations, rebuild TensorFlow with the appropriate compiler flags.
===>label 0
1/1 [==============================] - 0s 292ms/step
1/1 [==============================] - 0s 19ms/step
1/1 [==============================] - 0s 20ms/step
1/1 [==============================] - 0s 18ms/step
1/1 [==============================] - 0s 19ms/step
1/1 [==============================] - 0s 19ms/step
1/1 [==============================] - 0s 19ms/step
1/1 [==============================] - 0s 19ms/step
1/1 [==============================] - 0s 21ms/step
1/1 [==============================] - 0s 19ms/step
1/1 [==============================] - 0s 23ms/step
...
```

# Training

We use 68k PBMCs scRNA-seq data downloaded from [https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis](https://github.com/10XGenomics/single-cell-3prime-paper/tree/master/pbmc68k_analysis). It has the content.

* [x] `fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar`
* [x] `68k_pbmc_barcodes_annotation.tsv` 



We use an in-house data processing method for filtering cells and genes that are not proper for analysis. This method can be programmatically accessed with `uc.mat.data_in`. We can use it to read the data and further get the portions for training and testing.

:material-language-python: `Python`
``` py linenums="1"
import umiche as uc

x_train, y_train, x_test, y_test = uc.mat.data_in(
    params={
        'data_fpn': to('data/deepconvcvae/flt.h5'),
        'data_rv_zero_fpn': to('data/deepconvcvae/flt_rv_zero_col.h5'),
    },
    image_size=144,
)
```

:material-console: `console`
``` shell
Train: (41147,) Test: (27432,)
Train: (41147,) Test: (27432,)
Train: (41147,) Test: (27432,)
Train: (41147,) Test: (27432,)
Train: (41147,) Test: (27432,)
['CD8+ Cytotoxic T', 'CD8+/CD45RA+ Naive Cytotoxic', 'CD4+/CD45RO+ Memory', 'CD19+ B', 'CD4+/CD25 T Reg', 'CD56+ NK', 'CD4+ T Helper2', 'CD4+/CD45RA+/CD25- Naive T', 'CD34+', 'Dendritic', 'CD14+ Monocyte']
(41147, 144, 144, 1)
(41147,)
(27432, 144, 144, 1)
(27432,)
```

The following code will allow us to train the ConvCVAE framework. We can actually use any other types of training data, provided that you have tailored it to suit the input format (i.e. `numpy.ndarray`).

:material-language-python: `Python`
``` py linenums="1"
import umiche as uc

uc.mat.train_deepconvcvae(
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
)
```

:material-console: `console`
``` shell
KerasTensor(type_spec=TensorSpec(shape=(None, 144, 144, 1), dtype=tf.float32, name='encoder_input'), name='encoder_input', description="created by layer 'encoder_input'")
KerasTensor(type_spec=TensorSpec(shape=(None, 11), dtype=tf.float32, name='class_labels'), name='class_labels', description="created by layer 'class_labels'")
2024-07-30 23:13:52.493677: I tensorflow/core/platform/cpu_feature_guard.cc:182] This TensorFlow binary is optimized to use available CPU instructions in performance-critical operations.
To enable the following instructions: SSE SSE2 SSE3 SSE4.1 SSE4.2 AVX AVX2 FMA, in other operations, rebuild TensorFlow with the appropriate compiler flags.
Model: "ccvae"
__________________________________________________________________________________________________
 Layer (type)                Output Shape                 Param #   Connected to                  
==================================================================================================
 encoder_input (InputLayer)  [(None, 144, 144, 1)]        0         []                            
                                                                                                  
 class_labels (InputLayer)   [(None, 11)]                 0         []                            
                                                                                                  
 encoder (Functional)        [(None, 2),                  917412    ['encoder_input[0][0]',       
                              (None, 2),                             'class_labels[0][0]']        
                              (None, 2)]                                                          
                                                                                                  
 decoder (Functional)        (None, 144, 144, 1)          586465    ['encoder[0][2]',             
                                                                     'class_labels[0][0]']        
                                                                                                  
 dense (Dense)               (None, 20736)                248832    ['class_labels[0][0]']        
                                                                                                  
 reshape (Reshape)           (None, 144, 144, 1)          0         ['dense[0][0]']               
                                                                                                  
 concatenate (Concatenate)   (None, 144, 144, 2)          0         ['encoder_input[0][0]',       
                                                                     'reshape[0][0]']             
                                                                                                  
 conv2d (Conv2D)             (None, 72, 72, 16)           304       ['concatenate[0][0]']         
                                                                                                  
 conv2d_1 (Conv2D)           (None, 36, 36, 32)           4640      ['conv2d[0][0]']              
                                                                                                  
 flatten (Flatten)           (None, 41472)                0         ['conv2d_1[0][0]']            
                                                                                                  
 dense_1 (Dense)             (None, 16)                   663568    ['flatten[0][0]']             
                                                                                                  
 z_log_var (Dense)           (None, 2)                    34        ['dense_1[0][0]']             
                                                                                                  
 z_mean (Dense)              (None, 2)                    34        ['dense_1[0][0]']             
                                                                                                  
 tf.reshape_1 (TFOpLambda)   (None,)                      0         ['decoder[0][0]']             
                                                                                                  
 tf.reshape (TFOpLambda)     (None,)                      0         ['encoder_input[0][0]']       
                                                                                                  
 tf.__operators__.add (TFOp  (None, 2)                    0         ['z_log_var[0][0]']           
 Lambda)                                                                                          
                                                                                                  
 tf.math.square (TFOpLambda  (None, 2)                    0         ['z_mean[0][0]']              
 )                                                                                                
                                                                                                  
 tf.convert_to_tensor (TFOp  (None,)                      0         ['tf.reshape_1[0][0]']        
 Lambda)                                                                                          
                                                                                                  
 tf.cast (TFOpLambda)        (None,)                      0         ['tf.reshape[0][0]']          
                                                                                                  
 tf.math.subtract (TFOpLamb  (None, 2)                    0         ['tf.__operators__.add[0][0]',
 da)                                                                 'tf.math.square[0][0]']      
                                                                                                  
 tf.math.exp (TFOpLambda)    (None, 2)                    0         ['z_log_var[0][0]']           
                                                                                                  
 tf.math.squared_difference  (None,)                      0         ['tf.convert_to_tensor[0][0]',
  (TFOpLambda)                                                       'tf.cast[0][0]']             
                                                                                                  
 tf.math.subtract_1 (TFOpLa  (None, 2)                    0         ['tf.math.subtract[0][0]',    
 mbda)                                                               'tf.math.exp[0][0]']         
                                                                                                  
 tf.math.reduce_mean (TFOpL  ()                           0         ['tf.math.squared_difference[0
 ambda)                                                             ][0]']                        
                                                                                                  
 tf.math.reduce_sum (TFOpLa  (None,)                      0         ['tf.math.subtract_1[0][0]']  
 mbda)                                                                                            
                                                                                                  
 tf.math.multiply (TFOpLamb  ()                           0         ['tf.math.reduce_mean[0][0]'] 
 da)                                                                                              
                                                                                                  
 tf.math.multiply_1 (TFOpLa  (None,)                      0         ['tf.math.reduce_sum[0][0]']  
 mbda)                                                                                            
                                                                                                  
 tf.__operators__.add_1 (TF  (None,)                      0         ['tf.math.multiply[0][0]',    
 OpLambda)                                                           'tf.math.multiply_1[0][0]']  
                                                                                                  
 tf.math.reduce_mean_1 (TFO  ()                           0         ['tf.__operators__.add_1[0][0]
 pLambda)                                                           ']                            
                                                                                                  
 add_loss (AddLoss)          ()                           0         ['tf.math.reduce_mean_1[0][0]'
                                                                    ]                             
                                                                                                  
==================================================================================================
Total params: 1503877 (5.74 MB)
Trainable params: 1503877 (5.74 MB)
Non-trainable params: 0 (0.00 Byte)
__________________________________________________________________________________________________
  18/1286 [..............................] - ETA: 2:58 - loss: 3376.0042
```