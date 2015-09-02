
import sys

headers = sys.stdin.readline().strip().split('\t')
for idx, line in enumerate(sys.stdin):
  print "===== Sample %i =====" % idx
  fields = line.strip('\n').split('\t')
  for jdx, field in enumerate(fields):
    print "%24s: %s" % ( headers[jdx], field )
