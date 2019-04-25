#############################################################################################
# Scores when setting the prediction of the most likely outcome to 1 according to the pdf's, with one-sided pdfs #
#############################################################################################

# iter = 150, bin_width = 0.1, factor = 1
# 0.5744 (with most likely)

# 0.4796 (with probabilities)

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
factor = 1
bins = 100
bin_width = 0.1

TRAIN_ON_EVERYTHING = True

from pandas import *
from numpy import *
import scipy.optimize
from random import *
import datetime
from matplotlib.pyplot import *
import sys


#############
# Functions #
#############

def getrank(team):
    return rank.ix[index[team]].ix['Rank']

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
    if not 'Sign' in data:
        data.insert(len(data.columns), 'Sign', zeros(len(data)))
    for match in data.index:
        t1 = data.get_value(match,'Team1')
        t2 = data.get_value(match,'Team2')
        rank1 = rank.get_value(index[t1],'Rank')
        rank2 = rank.get_value(index[t2],'Rank')
        rankdiff = abs(rank1-rank2)
        data.RankDiff.ix[match] = rankdiff
        data.Sign.ix[match] = sign(rank1-rank2)

    #Criteria for losing/winning/tieing with no hometeams involved (being the highest ranked)
    wins = data[((data.Sign == 1) & (data.Outcome == 1.0) & (data.HomeTeam == 0)) | ((data.Sign == -1) & (data.Outcome == 0.0) & (data.HomeTeam == 0))]
    ties = data[(data.Outcome == 0.5) & (data.HomeTeam == 0)]
    losses = data[((data.Sign == 1) & (data.Outcome == 0.0) & (data.HomeTeam == 0)) | ((data.Sign == -1) & (data.Outcome == 1.0) & (data.HomeTeam == 0))]


    # MUST BE A POWER OF TWO
    binlist = [bin_width*i for i in range(bins)]
    #binlist = [x - binlist[bins/2] for x in binlist]
    #print binlist
    # Dictionary from rank diff to bin index
    global bin_index
    bin_index = dict(zip(binlist, range(len(binlist))))
    w = hist(wins.RankDiff, bins=binlist,   normed=False)
    t = hist(ties.RankDiff, bins=binlist,   normed=False)
    l = hist(losses.RankDiff, bins=binlist, normed=False)

    #Criteria for losing/winning/tieing as a hometeam (being the highest ranked)
    wins = data[((data.Sign == 1) & (data.Outcome == 1.0) & (data.HomeTeam == 1))]
    ties = data[((data.Sign == 1) & (data.Outcome == 0.5) & (data.HomeTeam == 1))]
    losses = data[((data.Sign == 1) & (data.Outcome == 0.0) & (data.HomeTeam == 1))]

    wh = hist(wins.RankDiff, bins=binlist,   normed=False)
    th = hist(ties.RankDiff, bins=binlist,   normed=False)
    lh = hist(losses.RankDiff, bins=binlist, normed=False)

    #Criteria for losing/winning/tieing against a hometeam (being the highest ranked)
    wins = data[((data.Sign == -1) & (data.Outcome == 0.0) & (data.HomeTeam == 1))]
    ties = data[((data.Sign == -1) & (data.Outcome == 0.5) & (data.HomeTeam == 1))]
    losses = data[((data.Sign == -1) & (data.Outcome == 1.0) & (data.HomeTeam == 1))]


    wah = hist(wins.RankDiff, bins=binlist,   normed=False)
    tah = hist(ties.RankDiff, bins=binlist,   normed=False)
    lah = hist(losses.RankDiff, bins=binlist, normed=False)


    return (w,t,l,wh,th,lh,wah,tah,lah)

#Prediction algorithm, take most likely outcome according to rankdiff and set prediction to 1
def prediction(data,rank,pdf_nohome,pdf_home,pdf_ahome):

    resultcsv = []
    pred_value = 0.0
    bad_matches = []
    good_matches = []

    (nohome, home, ahome) = map(fit_pdf, [pdf_nohome, pdf_home, pdf_ahome])
    for match in data.index:
        t1 = data.get_value(match,'Team1')
        t2 = data.get_value(match,'Team2')
        rank1 = rank.get_value(index[t1],'Rank')
        rank2 = rank.get_value(index[t2],'Rank')
        rankdiff = abs(rank1-rank2)
        signn = sign(rank1-rank2)

        #Find the right bin
        binn = floor(rankdiff/bin_width)*bin_width
        binn = bin_index[binn]

        #Choose pdf
        if signn == 1 and data.get_value(match,'Outcome') == 1 and data.get_value(match,'HomeTeam') == 1:
            (a,b,c) = home(rankdiff)
        elif signn == 1 and data.get_value(match,'Outcome') == 0.5 and data.get_value(match,'HomeTeam') == 1:
            (a,b,c) = home(rankdiff)
        elif signn == 1 and data.get_value(match,'Outcome') == 0 and data.get_value(match,'HomeTeam') == 1:
            (a,b,c) = home(rankdiff)
        elif signn == -1 and data.get_value(match,'Outcome') == 1 and data.get_value(match,'HomeTeam') == 1:
            (a,b,c) = ahome(rankdiff)
        elif signn == -1 and data.get_value(match,'Outcome') == 0.5 and data.get_value(match,'HomeTeam') == 1:
            (a,b,c) = ahome(rankdiff)
        elif signn == -1 and data.get_value(match,'Outcome') == 0 and data.get_value(match,'HomeTeam') == 1:
            (a,b,c) = ahome(rankdiff)
        else:
            (a,b,c) = nohome(rankdiff)

        #Predict outcome as most likely bin (convert from highest rank perspective to team1 perspective!)
        if signn == 1:
            predicted_outcome = 1 - argmax([a,b,c]) / 2.0
        else:
            predicted_outcome = 1 - (1 - argmax([a,b,c]) / 2.0)

        #If match, add (This is for self scoring)
        if abs(data.get_value(match,'Outcome') - predicted_outcome) < 0.01:
            pred_value += 1.0

        resultcsv.append([0,0,0])
        resultcsv[-1][int(1-predicted_outcome) * 2] = 1

    # with open('result.csv', 'w') as f:
        # for w,t,l in resultcsv:
            # f.write("%s,%s,%s\n" % (w,t,l))
        # #Predict outcome with probabilities
        # tmp = [a,b,c]
        # if signn == 1:
        #     pred_value += tmp[int((1-data.get_value(match,'Outcome'))*2)]
        # else:
        #     pred_value += tmp[int(data.get_value(match,'Outcome')*2)]


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
        # resultcsv.append(pdf[binn][1:])
    with open('result.csv', 'w') as f:
        for w,t,l in resultcsv:
            f.write("%s,%s,%s\n" % (w,t,l))
    return pred_value/len(data)


def fit_pdf(pdf):
    bins = pdf[:,0]
    wins = pdf[:,1]
    ties = pdf[:,2]
    loss = pdf[:,3]
    wins_model = lambda x, a, b, c: 1.0 / (a + b * exp(c*x))
    ties_model = lambda x, b, c: b * exp(c*x)
    loss_model = lambda x, b, c: b * exp(c*x)
    (w, cov) = scipy.optimize.curve_fit(wins_model, bins, wins)
    (t, cov) = scipy.optimize.curve_fit(ties_model, bins, ties)
    (l, cov) = scipy.optimize.curve_fit(loss_model, bins, loss)
    wins_model = lambda x: 1.0 / (w[0] + w[1] * exp(w[2]*x))
    ties_model = lambda x: t[0] * exp(t[1]*x)
    loss_model = lambda x: l[0] * exp(l[1]*x)
    hb = bin_width * 0.5
    def applymodel(x):
        vals = array((wins_model(x+hb), ties_model(x+hb), loss_model(x+hb)))
        return tuple(vals / sum(vals))
    return (lambda x: applymodel(x))


def home(x):
    return (1/(0.993745588217890352 +
 1.20447136550852012*exp(-1.16527348239199668*(x-0.5*bin_width)))/(1/(0.993745588217890352 +
 1.20447136550852012*exp(-1.16527348239199668*(x-0.5*bin_width)))+
0.411372791788055856*exp(-0.827794997859682452*(x-0.5*bin_width))+
0.179762416000732405*exp(-0.975266212675918470*(x-0.5*bin_width))),
0.411372791788055856*exp(-0.827794997859682452*(x-0.5*bin_width))/(1/(0.993745588217890352 +
 1.20447136550852012*exp(-1.16527348239199668*(x-0.5*bin_width)))+
0.411372791788055856*exp(-0.827794997859682452*(x-0.5*bin_width))+
0.179762416000732405*exp(-0.975266212675918470*(x-0.5*bin_width))),
0.179762416000732405*exp(-0.975266212675918470*(x-0.5*bin_width))/(1/(0.993745588217890352 +
 1.20447136550852012*exp(-1.16527348239199668*(x-0.5*bin_width)))+
0.411372791788055856*exp(-0.827794997859682452*(x-0.5*bin_width))+
0.179762416000732405*exp(-0.975266212675918470*(x-0.5*bin_width))))


def nohome(x):
    return (1/(0.999381312615316369 +
 1.60047731001837545*exp(-1.06709836157930476*(x-0.5*bin_width)))/
(1/(0.999381312615316369 +
 1.60047731001837545*exp(-1.06709836157930476*(x-0.5*bin_width)))+
0.363029702153055662*exp(-0.632174997542004040 *(x-0.5*bin_width))+
0.326256890083060776*exp(-0.906662594450005946 *(x-0.5*bin_width))),
0.363029702153055662*exp(-0.632174997542004040 *(x-0.5*bin_width))/
(1/(0.999381312615316369 +
 1.60047731001837545*exp(-1.06709836157930476*(x-0.5*bin_width)))+
0.363029702153055662*exp(-0.632174997542004040 *(x-0.5*bin_width))+
0.326256890083060776*exp(-0.906662594450005946 *(x-0.5*bin_width))),
0.326256890083060776*exp(-0.906662594450005946 *(x-0.5*bin_width))/
(1/(0.999381312615316369 +
 1.60047731001837545*exp(-1.06709836157930476*(x-0.5*bin_width)))+
0.363029702153055662*exp(-0.632174997542004040 *(x-0.5*bin_width))+
0.326256890083060776*exp(-0.906662594450005946 *(x-0.5*bin_width))))


def ahome1(x):
    return (1/(0.974568627036329208 +
 2.65013610774447005*exp(-0.835449469892954494 *(x-0.5*bin_width)))/
 (1/(0.974568627036329208 +
 2.65013610774447005*exp(-0.835449469892954494 *(x-0.5*bin_width)))+
0.398367959209470331*exp(-0.458448718924871196 *(x-0.5*bin_width))+
0.456455748206617582*exp(-0.691275383705087610 *(x-0.5*bin_width))),
0.398367959209470331*exp(-0.458448718924871196 *(x-0.5*bin_width))/
 (1/(0.974568627036329208 +
 2.65013610774447005*exp(-0.835449469892954494 *(x-0.5*bin_width)))+
0.398367959209470331*exp(-0.458448718924871196 *(x-0.5*bin_width))+
0.456455748206617582*exp(-0.691275383705087610 *(x-0.5*bin_width))),
0.456455748206617582*exp(-0.691275383705087610 *(x-0.5*bin_width))/
 (1/(0.974568627036329208 +
 2.65013610774447005*exp(-0.835449469892954494 *(x-0.5*bin_width)))+
0.398367959209470331*exp(-0.458448718924871196 *(x-0.5*bin_width))+
0.456455748206617582*exp(-0.691275383705087610 *(x-0.5*bin_width))))

#############################
# Initialisation and checks #
#############################

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

    if TRAIN_ON_EVERYTHING:
        data_training = data
        data_testing  = unplayed
    else:
        data_training = data[start:end].reset_index()
        data_testing  = data[end:moreend].reset_index()

    #Main loop
    # TODO Look at whether the distribution of ranks change as we iterate...
    rank_distribution = DataFrame(zeros((iter, len(rank.Team))), columns=rank.Team)
    for i in range(iter):
        print 'Iteration number: ', i

        #Compute neighbor averages
        # neighbor_avg(i,data,rank,index)

        #Learning rate
        eta = pow((1+0.1*iter)/(i+0.1*iter),0.602)
        # eta = 1

        #Create shuffled list
        x = range(len(data_training))
        shuffle(x)

        #Update step
        update_ranks(eta,data_training,rank,index,x)
        rank_distribution.ix[i] = rank.Rank
        # print rank[:10]

    #Create histograms
    foo = histograms(data_training,rank)

    #Pdf without distinction between hometeams or not (#NOT CORRECT!) (not being used)
    pdf = np.append(np.append(np.append(np.transpose([foo[0][1]]),np.transpose([np.append(map(float,foo[0][0]+foo[3][0])/(foo[0][0]+foo[1][0]+foo[2][0]+foo[3][0]+foo[4][0]+foo[5][0]),0)]),1),np.transpose([np.append(map(float,foo[1][0]+foo[4][0])/(foo[0][0]+foo[1][0]+foo[2][0]+foo[3][0]+foo[4][0]+foo[5][0]),0)]),1),np.transpose([np.append(map(float,foo[2][0]+foo[5][0])/(foo[0][0]+foo[1][0]+foo[2][0]+foo[3][0]+foo[4][0]+foo[5][0]),0)]),1)

    #Pdf for non-home team matches
    pdf_nohome = np.append(np.append(np.append(np.transpose([foo[0][1]]),np.transpose([np.append(map(float,foo[0][0])/(foo[0][0]+foo[1][0]+foo[2][0]),0)]),1),np.transpose([np.append(map(float,foo[1][0])/(foo[0][0]+foo[1][0]+foo[2][0]),0)]),1),np.transpose([np.append(map(float,foo[2][0])/(foo[0][0]+foo[1][0]+foo[2][0]),0)]),1)
    pdf_home   = np.append(np.append(np.append(np.transpose([foo[0][1]]),np.transpose([np.append(map(float,foo[3][0])/(foo[3][0]+foo[4][0]+foo[5][0]),0)]),1),np.transpose([np.append(map(float,foo[4][0])/(foo[3][0]+foo[4][0]+foo[5][0]),0)]),1),np.transpose([np.append(map(float,foo[5][0])/(foo[3][0]+foo[4][0]+foo[5][0]),0)]),1)
    pdf_ahome   = np.append(np.append(np.append(np.transpose([foo[0][1]]),np.transpose([np.append(map(float,foo[6][0])/(foo[6][0]+foo[7][0]+foo[8][0]),0)]),1),np.transpose([np.append(map(float,foo[7][0])/(foo[6][0]+foo[7][0]+foo[8][0]),0)]),1),np.transpose([np.append(map(float,foo[8][0])/(foo[6][0]+foo[7][0]+foo[8][0]),0)]),1)

    #Replace extreme rank differences with 1,0,0 or 0,0,1
    for i,r in enumerate(pdf):
        if any(isnan(r)):
            pdf[i][1:] = array([1,0,0])
    pdf = pdf[:-1]
    for i,r in enumerate(pdf_nohome):
        if any(isnan(r)):
            pdf_nohome[i][1:] = array([1,0,0])
    pdf_nohome = pdf_nohome[:-1]
    for i,r in enumerate(pdf_home):
        if any(isnan(r)):
            pdf_home[i][1:] = array([1,0,0])
    pdf_home = pdf_home[:-1]
    for i,r in enumerate(pdf_ahome):
        if any(isnan(r)):
            pdf_ahome[i][1:] = array([1,0,0])
    pdf_ahome = pdf_ahome[:-1]


    #Predict matches
    a = prediction(data_testing,rank,pdf_nohome,pdf_home,pdf_ahome)

    print a

def hometeam_advantage(teamname):
    home_ad_wins = Series([x and y and z for (x, y, z) in zip(data.Team1 == teamname,
                                                    data.HomeTeam == 1,
                                                    data.Score1 > data.Score2)])
    home_ad_ties = Series([x and y and z for (x, y, z) in zip(data.Team1 == teamname,
                                                    data.HomeTeam == 1,
                                                    data.Score1 == data.Score2)])
    home_ad_loss = Series([x and y and z for (x, y, z) in zip(data.Team1 == teamname,
                                                    data.HomeTeam == 1,
                                                    data.Score1 < data.Score2)])

    not_home_wins = Series([x and y and z for (x, y, z) in zip(data.Team1 == teamname,
                                              data.HomeTeam == 0,
                                              data.Score1 > data.Score2)])
    not_home_ties = Series([x and y and z for (x, y, z) in zip(data.Team1 == teamname,
                                              data.HomeTeam == 0,
                                              data.Score1 == data.Score2)])
    not_home_loss = Series([x and y and z for (x, y, z) in zip(data.Team1 == teamname,
                                              data.HomeTeam == 0,
                                              data.Score1 < data.Score2)])

    (w,t,l) = len(data[home_ad_wins]), len(data[home_ad_ties]), len(data[home_ad_loss])
    print w / float(w+t+l), t / float(w+t+l), l / float(w+t+l)
    (w,t,l) = len(data[not_home_wins]), len(data[not_home_ties]), len(data[not_home_loss])
    print w / float(w+t+l), t / float(w+t+l), l / float(w+t+l)
    return home_ad_wins



