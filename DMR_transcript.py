//calculate the distance of differentially methylated site (DMR) and transcript in reference genome

#!/usr/bin/python3
import pandas
import matplotlib.pyplot as plt
 
df = pandas.read_table("/Users/lucydiamondsky/genomic/v27lift37.gpd",header=None)
of = pandas.read_table("/Users/lucydiamondsky/genomic/DMR_new2.txt",header=None)
df.ix[:,:1].groupby(1).count().plot.bar(legend=False)
plt.show()

result = {}
for x in range(0,194):
    dmr =  (of.ix[x,1] + of.ix[x,2])/2
    for i in range(0,202697):
        distance = 1000000
        if of.ix[x,0] == df.ix[i,1]:
            distance2 = abs(dmr - df.ix[i,3])
            distance3 = abs(dmr - df.ix[i,4])
            if distance2 < distance:
                distance = distance2
                if distance3 <distance:
                    distance = distance3
                    result[dmr] = [distance,df.ix[i,0]]

result = {}
for x in range(0,5):
    dmr =  (of.ix[x,1] + of.ix[x,2])/2
    for i in range(0,202696):
        distance = 1000000
        if of.ix[x,0] == df.ix[i,1]:
            if dmr > df.ix[i,3] and dmr < df.ix[i,4]:
                result[dmr] = ['the dmr is within transcript',df.ix[i,0]
            else:
                distance2 = abs(dmr - df.ix[i,3])
                distance3 = abs(dmr - df.ix[i,4])
                if distance2 < distance:
                    distance = distance2
                    if distance3 <distance:
                        distance = distance3
                        result[dmr] = [distance,df.ix[i,0]]