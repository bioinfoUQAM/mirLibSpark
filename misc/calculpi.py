from pyspark import SparkConf, SparkContext
import random

def inside(p):
    x, y = random.random(), random.random()
    return x*x + y*y < 1


conf = SparkConf()
sc = SparkContext(conf = conf)
NUM_SAMPLES=10000

count = sc.parallelize(xrange(0, NUM_SAMPLES)) \
             .filter(inside).count()

print "Pi is roughly %f" % (4.0 * count / NUM_SAMPLES)