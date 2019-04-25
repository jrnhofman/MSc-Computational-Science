from pandas import *
from numpy import *
from random import *
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
iter = 60
smoothness = 1
factor = 20

#############
# Functions #
#############

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


def pred_outcome(Team1, Team2, rank):
    rankA = rank.get_value(index[Team1],'Rank')
    rankB = rank.get_value(index[Team2],'Rank')
    sgn = sign(rankB - rankA)
    # return 1/(1 + exp(smoothness * sgn * sqrt(abs(rankB - (rankA + gamma)))))
    return 1/(1 + exp(smoothness * (rankB - (rankA + gamma))))

def spread(rank):
    ret = []
    rank5 = rank.sort('Rank', ascending=True)
    rank5 = rank5.reset_index()
    test_teams = rank5.Team
    for team1,team2 in zip(test_teams, test_teams[1:]):
        ret.append(pred_outcome(team1, team2, rank))
    return ret


#############################
# Initialisation and checks #
#############################

#Read in file
data = read_csv('FIFA0212.csv',parse_dates=True)

#Throw away unplayed matches and unused attributes
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
# print rank[:10]


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
    wins = data[(data.Outcome == 1) & (data.HomeTeam == 0)]
    ties = data[data.Outcome == 0.5]
    losses = data[data.Outcome == 0.0]
    # w = hist((wins.RankDiff, ties.RankDiff, losses.RankDiff), bins=50, normed=True)
    w = hist(wins.RankDiff, bins=50,   normed=False)
    t = hist(ties.RankDiff, bins=50,   normed=False)
    l = hist(losses.RankDiff, bins=50, normed=False)
    return (w,t,l)

#Main loop
segments = 10
totals = 0
for piece in range(segments-1):
    start = int(len(data) * piece / 10.)
    end   = int(len(data) * (piece+1) / 10.)
    moreend = int(len(data) * (piece+2) / 10.)

    data_training = data[start:end].reset_index()
    data_testing  = data[end:moreend].reset_index()

    for i in range(iter):
        # print 'Iteration number: ', i

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

    total = 0
    successful_tests = 0
    total_tests = len(data_testing)
    for match in data_testing.index:
        t1 = data_testing.get_value(match,'Team1')
        t2 = data_testing.get_value(match,'Team2')
        if not (any(data_training.Team1 == t1) or
                any(data_training.Team2 == t1) and
                any(data_training.Team1 == t2) or
                any(data_training.Team2 == t2)):
                # print data_testing.ix[match]
                continue
        successful_tests += 1
        prediction = round(pred_outcome(t1, t2, rank)*2.0 ) / 2.0
        outcome = data_testing.get_value(match,'Outcome')
        rank1 = rank.get_value(index[t1],'Rank')
        rank2 = rank.get_value(index[t2],'Rank')
        print "Team1: (%s, %f)  Team2: (%s, %f)\nOutcome: %s, Prediction:  %s" % (t1,
                                                                        rank1,
                                                                        t2,
                                                                        rank2,
                                                                        outcome,
                                                                        prediction)
        if raw_input() == 'q\n':
            sys.exit(0)
        total += abs(prediction - outcome)
    total /= successful_tests
    totals += totals
    print piece, successful_tests, total
totals /= segments
print totals

# print rank[:25]