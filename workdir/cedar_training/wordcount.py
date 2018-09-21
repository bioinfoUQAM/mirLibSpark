from __future__ import print_function

import sys
from operator import add

from pyspark import SparkContext


partition = 4 #= set as the same number as the processors


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: wordcount <file>", file=sys.stderr)
        exit(-1)
    sc = SparkContext(appName="PythonWordCount")
    lines = sc.textFile(sys.argv[1], partition)
    counts = lines.flatMap(lambda x: x.split(' ')) \
                  .map(lambda x: (x, 1)) \
                  .reduceByKey(add)
    output = counts.collect()
    for (word, count) in output:
        print("%s: %i" % (word, count))

    sc.stop()

#= https://github.com/eBay/Spark/blob/master/examples/src/main/python/wordcount.py
# spark-submit wordcount.py newyorktimes.txt
