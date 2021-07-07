print(sys.version)

import os
print(os.environ['CUDA_VISIBLE_DEVICES'])

import tensorflow as tf
#print(tf.test.is_gpu_available())
print(tf.config.list_physical_devices('GPU'))

