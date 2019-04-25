#############################################################################################
# Scores when setting the prediction of the most likely outcome to 1 according to the pdf's #
#############################################################################################

# 0.5657 with home
# 0.52 without home

from pandas import *
from numpy import *
from random import *
import datetime
from matplotlib.pyplot import *
import sys

#########################
# Important assumptions #
#########################

# Actual scores don't matter, just win, loss or tie
# Home team is not incorporated yet, assumed that every Team1 has home advantage

##############
# Parameters #
##############

# gamma has no effect now.
gamma = 0.0
# lambd = 1e-7
lambd = 0.0
iter = 300
smoothness = 0.5
factor = 20
bins = 65
bin_width = 0.25

#############
# Functions #
#############

#Read dates
def readdate(s):
    return datetime.datetime.strptime(s, '%d/%m/%Y')

#Compute time weights
def compute_weights(data):

    #Weights given per month
    tmin = factor * 1./12 #Jan 2002
    tmax = factor * 10.4167 #May 2012

    data2 = zeros((len(data),2))

    for i in range(len(data)):
        string = data.get_value(i,'Date')[3:]
        sum = int((string[0:2]))*tmin + (int(string[3:])-2002) * tmin * 12 #based on month/year only
        data2[i][0] = ((1 + sum - tmin)/(1 + tmax - tmin))**2 #3.1 of the paper

        #Compute outcome scores
        if data.get_value(i,'Score1') > data.get_value(i,'Score2'):
            data2[i][1] = 1
        elif data.get_value(i,'Score1') < data.get_value(i,'Score2'):
            data2[i][1] = 0
        else:
            data2[i][1] = 0.5

    return data2

#Compute neighbor_averages according to algorithm from paper
def neighbor_avg(i,data,rank,ranks):

    #First iteration, compute sum of weights of played matches and neighborhood size
    if(i==0):
        nb_avg = zeros(len(rank))
        nb_size = zeros((len(rank),len(rank)))

        for j in range(len(data)):
            A = data.get_value(j,'Team1')
            B = data.get_value(j,'Team2')
            nb_avg[index[A]] += data.get_value(j,'Weight')
            nb_avg[index[B]] += data.get_value(j,'Weight')
            nb_size[index[A],index[B]] = 1
            nb_size[index[B],index[A]] = 1

        for j in range(len(rank)):
            rank.set_value(j,'Nb_avg_denom',nb_avg[j])
            rank.set_value(j,'Nb_size',nb_size.sum(axis=1)[j])

    #Every iteration, compute neighbor average
    nb_avg = zeros(len(rank))

    for j in range(len(data)):
        A = data.get_value(j,'Team1')
        B = data.get_value(j,'Team2')
        nb_avg[index[A]] += data.get_value(j,'Weight') * rank.get_value(index[B],'Rank')
        nb_avg[index[B]] += data.get_value(j,'Weight') * rank.get_value(index[A],'Rank')

    for j in range(len(rank)):
        rank.set_value(j,'Nb_avg',nb_avg[j]/rank.get_value(j,'Nb_avg_denom'))

#Ranking algorithm
def update_ranks(eta,data,rank,index,x):

    #For all shuffled matches, update rank
    for i in x:
        A = data.get_value(i,'Team1')
        B = data.get_value(i,'Team2')
        rankA = rank.get_value(index[A],'Rank')
        rankB = rank.get_value(index[B],'Rank')

        pred_outcome = 1/(1 + exp(rankB - (rankA + gamma)))
        outcome = data.get_value(i,'Outcome')

        rank.set_value(index[A],'Rank', rankA - eta *
                       (data.get_value(i,'Weight') * (pred_outcome - outcome)
                        * pred_outcome * (1 - pred_outcome)
                        # + (lambd/rank.get_value(index[A],'Nb_size'))
                            # * (rankA - rank.get_value(index[A],'Nb_avg'))
                       ))

        rank.set_value(index[B],'Rank', rankB - eta *
                       (-data.get_value(i,'Weight') * (pred_outcome - outcome)
                        * pred_outcome * (1 - pred_outcome)
                        # + (lambd/rank.get_value(index[B],'Nb_size'))
                            # * (rankB - rank.get_value(index[B],'Nb_avg'))
                       ))


#Predict outcome given teams
def pred_outcome(Team1, Team2, rank):
    rankA = rank.get_value(index[Team1],'Rank')
    rankB = rank.get_value(index[Team2],'Rank')
    sgn = sign(rankB - rankA)
    # return 1/(1 + exp(smoothness * sgn * sqrt(abs(rankB - (rankA + gamma)))))
    return 1/(1 + exp(smoothness * (rankB - (rankA + gamma))))

#Compute pred_outcome for neighbor teams in ranking list (should be around 0.5)
def spread(rank):
    ret = []
    rank5 = rank.sort('Rank', ascending=True)
    rank5 = rank5.reset_index()
    test_teams = rank5.Team
    for team1,team2 in zip(test_teams, test_teams[1:]):
        ret.append(pred_outcome(team1, team2, rank))
    return ret

#Create pdfs for home and non-home team matches
def histograms(data, rank):
    if not 'RankDiff' in data:
        data.insert(len(data.columns), 'RankDiff', zeros(len(data)))
    for match in data.index:
        t1 = data.get_value(match,'Team1')
        t2 = data.get_value(match,'Team2')
        rank1 = rank.get_value(index[t1],'Rank')
        rank2 = rank.get_value(index[t2],'Rank')
        rankdiff = rank1-rank2
        data.RankDiff.ix[match] = rankdiff

    wins = data[(data.Outcome == 1.0) & (data.HomeTeam == 0)]
    ties = data[(data.Outcome == 0.5) & (data.HomeTeam == 0)]
    losses = data[(data.Outcome == 0.0) & (data.HomeTeam == 0)]

    # MUST BE A POWER OF TWO
    binlist = [bin_width*i for i in range(bins)]
    binlist = [x - binlist[bins/2] for x in binlist]
    # Dictionary from rank diff to bin index
    global bin_index
    bin_index = dict(zip(binlist, range(len(binlist))))
    w = hist(wins.RankDiff, bins=binlist,   normed=False)
    t = hist(ties.RankDiff, bins=binlist,   normed=False)
    l = hist(losses.RankDiff, bins=binlist, normed=False)

    wins = data[(data.Outcome == 1.0) & (data.HomeTeam == 1)]
    ties = data[(data.Outcome == 0.5) & (data.HomeTeam == 1)]
    losses = data[(data.Outcome == 0.0) & (data.HomeTeam == 1)]

    wh = hist(wins.RankDiff, bins=binlist,   normed=False)
    th = hist(ties.RankDiff, bins=binlist,   normed=False)
    lh = hist(losses.RankDiff, bins=binlist, normed=False)

    return (w,t,l,wh,th,lh)

#Prediction algorithm, take most likely outcome according to rankdiff and set prediction to 1
def prediction(data,rank,pdf_nohome,pdf_home):

    resultcsv = []
    pred_value = 0.0
    bad_matches = []
    good_matches = []

    for match in data.index:
        t1 = data.get_value(match,'Team1')
        t2 = data.get_value(match,'Team2')
        rank1 = rank.get_value(index[t1],'Rank')
        rank2 = rank.get_value(index[t2],'Rank')
        rankdiff = rank1-rank2

        #Find the right bin
        binn = floor(rankdiff/bin_width)*bin_width
        binn = bin_index[binn]

        #Choose pdf
        if data.get_value(match,'HomeTeam') == 1:
            pdf = pdf_home
        else:
            pdf = pdf_nohome

        #Predict outcome as most likely bin
        predicted_outcome = 1 - argmax(pdf[binn][1:]) / 2.0

        #If match, add
        if abs(data.get_value(match,'Outcome') - predicted_outcome) < 0.01:
            pred_value += 1.0
            #good_matches.append(True)
            #bad_matches.append(False)
            #else:
            # print data.get_value(match,'Outcome'), pdf[binn][1:]
            #good_matches.append(False)
            #bad_matches.append(True)
        # if data.get_value(match,'Outcome') == 1:
            # pred_value += pdf[binn][1]
        # elif data.get_value(match,'Outcome') == 0.5:
            # pred_value += pdf[binn][2]
        # else:
            # pred_value += pdf[binn][3]
        resultcsv.append(pdf[binn][1:])
    # with open('result.csv', 'w') as f:
        # for w,t,l in resultcsv:
            # f.write("%s,%s,%s\n" % (w,t,l))

    return pred_value/len(data)


#############################
# Initialisation and checks #
#############################

#Read in file
data = read_csv('FIFA0212.csv',parse_dates=True)

#Throw away unplayed matches and unused attributes
# unplayed = data[(data['Score1'] == '?')]
# unplayed.Date = map(readdate, unplayed.Date)
# unplayed = unplayed[unplayed.Date <= datetime.datetime(2012,06,05)]
# unplayed = unplayed.join(DataFrame(range(len(unplayed)),columns=['Outcome']))
# data2 = compute_weights(unplayed)
# unplayed = unplayed.join(DataFrame(data2,columns=['Weight','Outcome']))

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


#############################
# The algorithm starts here #
#############################

#Add columns to ranking list for neighbor averages and neighborhood size (denom=denominator)
foo = zeros((len(rank),3))
rank = rank.join(DataFrame(foo,columns=['Nb_avg','Nb_avg_denom','Nb_size']))

#Compute weights and outcome score (1=win, 0=loss, 1/2=tie)
data2 = compute_weights(data)
data = data.join(DataFrame(data2,columns=['Weight','Outcome']))
# print data[:10].to_string()
# print data[-10:].to_string()

#Dictionary to compute index given team name
index = dict(zip(rank.Team,[i for i in range(len(rank))]))
bin_index = None
# print rank[:10]



segments = 2.
totals = 0

#Cross validation
for piece in range(int(segments)-1):
    start = int(len(data) * piece / segments)
    end   = int(len(data) * (piece+1) / segments)
    moreend = int(len(data) * (piece+2) / segments)

    data_training = data[start:end].reset_index()
    # data_training = data
    data_testing  = data[end:moreend].reset_index()
    # data_testing  = unplayed

    #Main loop
    for i in range(iter):
        print 'Iteration number: ', i

        #Compute neighbor averages
        # neighbor_avg(i,data,rank,index)

        #Learning rate
        eta = pow((1+0.1*iter)/(i+0.1*iter),0.602)

        #Create shuffled list
        x = range(len(data_training))
        shuffle(x)

        #Update step
        update_ranks(eta,data_training,rank,index,x)
        # print rank[:10]

    #Create histograms
    foo = histograms(data_training,rank)

    #Pdf for home team matches
    pdf = np.append(np.append(np.append(np.transpose([foo[0][1]]),np.transpose([np.append(map(float,foo[0][0]+foo[3][0])/(foo[0][0]+foo[1][0]+foo[2][0]+foo[3][0]+foo[4][0]+foo[5][0]),0)]),1),np.transpose([np.append(map(float,foo[1][0]+foo[4][0])/(foo[0][0]+foo[1][0]+foo[2][0]+foo[3][0]+foo[4][0]+foo[5][0]),0)]),1),np.transpose([np.append(map(float,foo[2][0]+foo[5][0])/(foo[0][0]+foo[1][0]+foo[2][0]+foo[3][0]+foo[4][0]+foo[5][0]),0)]),1)

    #Pdf for non-home team matches
    pdf_nohome = np.append(np.append(np.append(np.transpose([foo[0][1]]),np.transpose([np.append(map(float,foo[0][0])/(foo[0][0]+foo[1][0]+foo[2][0]),0)]),1),np.transpose([np.append(map(float,foo[1][0])/(foo[0][0]+foo[1][0]+foo[2][0]),0)]),1),np.transpose([np.append(map(float,foo[2][0])/(foo[0][0]+foo[1][0]+foo[2][0]),0)]),1)
    pdf_home   = np.append(np.append(np.append(np.transpose([foo[0][1]]),np.transpose([np.append(map(float,foo[3][0])/(foo[3][0]+foo[4][0]+foo[5][0]),0)]),1),np.transpose([np.append(map(float,foo[4][0])/(foo[3][0]+foo[4][0]+foo[5][0]),0)]),1),np.transpose([np.append(map(float,foo[5][0])/(foo[3][0]+foo[4][0]+foo[5][0]),0)]),1)

    #Replace extreme rank differences with 1,0,0 or 0,0,1
    for i,r in enumerate(pdf):
        if any(isnan(r)):
            if i < len(pdf)/2:
                pdf[i][1:] = array([0,0,1])
            else:
                pdf[i][1:] = array([1,0,0])
    pdf = pdf[:-1]
    for i,r in enumerate(pdf_nohome):
        if any(isnan(r)):
            if i < len(pdf_nohome)/2:
                pdf_nohome[i][1:] = array([0,0,1])
            else:
                pdf_nohome[i][1:] = array([1,0,0])
    pdf_nohome = pdf_nohome[:-1]
    for i,r in enumerate(pdf_home):
        if any(isnan(r)):
            if i < len(pdf_home)/2:
                pdf_home[i][1:] = array([0,0,1])
            else:
                pdf_home[i][1:] = array([1,0,0])
    pdf_home = pdf_home[:-1]

    #Predict matches
    a = prediction(data_testing,rank,pdf_nohome,pdf_home)
    print a


