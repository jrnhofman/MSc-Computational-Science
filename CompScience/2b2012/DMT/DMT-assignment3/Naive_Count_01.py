

from pandas import *
from numpy import *
from random import *
import datetime
from matplotlib.pyplot import *
import sys
TRAIN_ON_EVERYTHING = True


#Read dates
def readdate(s):
    return datetime.datetime.strptime(s, '%d/%m/%Y')

#Read in file
data = read_csv('FIFA0212.csv',parse_dates=True)

#Throw away unplayed matches and unused attributes
if TRAIN_ON_EVERYTHING:
    unplayed = data[(data['Score1'] == '?')]
    unplayed.Date = map(readdate, unplayed.Date)
    unplayed = unplayed[unplayed.Date <= datetime.datetime(2012,06,05)]
    unplayed = unplayed.join(DataFrame(range(len(unplayed)),columns=['Outcome']))

data = data[(data['Score1'] != '?')]
data = data[['Date','Team1','Team2','Score1','Score2','HomeTeam']]

#Throw away irrelevant countries
data = data[data['Team1'] != 'Yugoslavia']
data = data[data['Team2'] != 'Yugoslavia']
data = data[data['Team1'] != 'NetherlandsAntilles']
data = data[data['Team2'] != 'NetherlandsAntilles']
data = data[data['Team1'] != 'SerbiaandMontenegro']
data = data[data['Team2'] != 'SerbiaandMontenegro'].reset_index()

#Set datatypes
data['Score1'] = map(int,data['Score1'])
data['Score2'] = map(int,data['Score2'])
data['HomeTeam'] = map(int,data['HomeTeam'])

#Read in FIFA ranking list
rank = read_csv('Ranking.txt',sep='\t',header=None,names=['ActualRank','Team','Rank','bla2', 'bla3', 'bla4'] )
rank = rank[['Rank','Team']]
rank['Rank'] = map(float,rank['Rank'])
rank['Rank'] = zeros(len(rank))

#Checking whether spelling is the same in both lists
for i in range(len(rank)):
    for j in range(len(data)):
        if str(data.get_value(j,'Team1')) == rank.get_value(i,'Team') or data.get_value(j,'Team2') == rank.get_value(i,'Team'):
            break
        if j==len(data)-1:
            print 'Country not found!',rank.get_value(i,'Team')
# print 'All country names checked!'

data = data.join(DataFrame(zeros(len(data)),columns=['Outcome']))
data.Outcome[data.Score1 > data.Score2] = 1
data.Outcome[data.Score1 == data.Score2] = 0.5


resultcsv = []
for t1,t2 in zip(unplayed.Team1, unplayed.Team2):
    matches = data[((data.Team1 == t1) & (data.Team2 == t2))
       | ((data.Team1 == t2) & (data.Team2 == t1))]
    if len(matches) == 0:
        resultcsv.append([1/3.0]*3)
    else:
        w = matches[((matches.Outcome == 0) & (matches.Team1 == t2))
                    |((matches.Outcome == 1) & (matches.Team1 == t1))]
        t = matches[matches.Outcome == 0.5]
        w = len(w) / float(len(matches))
        t = len(t) / float(len(matches))
        l = 1 - w - t
        resultcsv.append([w,t,l])
with open('result.csv', 'w') as f:
    for w,t,l in resultcsv:
        f.write("%s,%s,%s\n" % (w,t,l))
