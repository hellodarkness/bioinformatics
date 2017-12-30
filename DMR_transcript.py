//calculate the distance of differentially methylated site (DMR) and transcript in reference genome
//http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=605193253_5uQOh6j9UyM6veEM8NZgpGxg7Lun
//transcript file achieved from reference genome
//browser session: http://epigenomegateway.wustl.edu/browser/?genome=hg19&session=EI1EV6BG4w&statusId=100728609
        
#!/usr/bin/python3
import pandas
import matplotlib.pyplot as plt

df = pandas.read_table(
    "/Users/lucydiamondsky/genomic/v27lift37.gpd", header=None)
of = pandas.read_table(
    "/Users/lucydiamondsky/genomic/DMR_new2.txt", header=None)
tf = pandas.read_table(
    "/Users/lucydiamondsky/genomic/tss", header=None)
df.ix[:, :1].groupby(1).count().plot.bar(legend=False)

#calculate distance to adjacent TSS
result3 = {}
for x in range(0, 194):
    dmr = (of.ix[x, 1] + of.ix[x, 2]) / 2
    distance = 1000000000
    for i in range(0, 68798):
        if (of.ix[x, 0] == tf.ix[i, 2]) and (tf.ix[i, 3] == '+'):
            distance2 = abs(of.ix[x, 1] - tf.ix[i, 4])
            distance3 = abs(of.ix[x, 2] - tf.ix[i, 4])
            distance4 = min(distance2, distance3)
            if distance4 < distance:
                distance = distance4
                result3[dmr] = [of.ix[x, 0], distance, tf.ix[i, 1]]
        if (of.ix[x, 0] == tf.ix[i, 2]) and (tf.ix[i, 3] == '-'):
            distance2 = abs(of.ix[x, 1] - tf.ix[i, 5])
            distance3 = abs(of.ix[x, 2] - tf.ix[i, 5])
            distance4 = min(distance2, distance3)
            if distance4 < distance:
                distance = distance4
                result3[dmr] = [of.ix[x, 0], distance, tf.ix[i, 1]]

for v in result3.values():
    print(format(v))

for v, k in result3.items():
    print('{v} {k}'.format(v=v, k=k))

#calculate distance to adjacent transcript
result = {}
for x in range(0, 195):
    dmr = (of.ix[x, 1] + of.ix[x, 2]) / 2
    distance = 1000000000
    for i in range(0, 202696):
        if of.ix[x, 0] == df.ix[i, 1]:
            if of.ix[x, 1] > df.ix[i, 3] and of.ix[x, 2] < df.ix[i, 4]:
                result[dmr] = [of.ix[x, 0],
                               'the dmr is located within transcript', df.ix[i, 0]]
            else:
                distance2 = abs(dmr - df.ix[i, 3])
                distance3 = abs(dmr - df.ix[i, 4])
                if distance2 < distance:
                    distance = distance2
                if distance3 < distance:
                    distance = distance3
                    result[dmr] = [of.ix[x, 0], distance, df.ix[i, 0]]

result1 = {}
for x in range(0, 195):
    dmr = (of.ix[x, 1] + of.ix[x, 2]) / 2
    for i in range(0, 202696):
        if of.ix[x, 0] == df.ix[i, 1]:
            if of.ix[x, 2] > df.ix[i, 3] and of.ix[x, 1] < df.ix[i, 3]:
                result1[dmr] = [of.ix[x, 0],
                                'the dmr is overlap with transcript', df.ix[i, 0]]
            if of.ix[x, 2] > df.ix[i, 4] and of.ix[x, 1] < df.ix[i, 4]:
                result1[dmr] = [of.ix[x, 0],
                                'the dmr is overlap with transcript', df.ix[i, 0]]

for v in result3.values():
    print(format(v))

for v, k in result3.items():
    print('{v} {k}'.format(v=v, k=k))

transcript = {}
for x in range(0, 194):
    result3.values()[x][2]
    transcript[x] = [result.values()[x][2]]
