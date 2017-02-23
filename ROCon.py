#!/usr/bin/env python

#based on plotroc.py

#!/usr/bin/env python
#This tool allow users to get SVM-prob ROC curve data for plotting by their preferred means
#Be sure to use properly scaled data

from svmutil import *
from sys import argv, platform
from os import path, popen
from random import randrange , seed
from operator import itemgetter
from time import sleep

def get_pos_deci(train_y, train_x, test_y, test_x, param):
        model = svm_train(train_y, train_x, param)
        #predict and grab decision value, assure deci>0 for label+,
        #the positive descision value = val[0]*labels[0]
        labels = model.get_labels()
        py, evals, deci = svm_predict(test_y, test_x, model)
        deci = [labels[0]*val[0] for val in deci]
        return deci,model

#input raw attributes, labels, param, cv_fold in decision value building
#output list of decision value, remember to seed(0)
def get_cv_deci(prob_y, prob_x, param, nr_fold):
        if nr_fold == 1 or nr_fold==0:
                deci,model = get_pos_deci(prob_y, prob_x, prob_y, prob_x, param)
                return deci
        deci, model = [], []
        prob_l = len(prob_y)

        #random permutation by swapping i and j instance
        for i in range(prob_l):
                j = randrange(i,prob_l)
                prob_x[i], prob_x[j] = prob_x[j], prob_x[i]
                prob_y[i], prob_y[j] = prob_y[j], prob_y[i]

        #cross training : folding
        for i in range(nr_fold):
                begin = i * prob_l // nr_fold
                end = (i + 1) * prob_l // nr_fold
                train_x = prob_x[:begin] + prob_x[end:]
                train_y = prob_y[:begin] + prob_y[end:]
                test_x = prob_x[begin:end]
                test_y = prob_y[begin:end]
                subdeci, submdel = get_pos_deci(train_y, train_x, test_y, test_x, param)
                deci += subdeci
        return deci

#processing argv and set some global variables
def proc_argv(argv = argv):
        #print("Usage: %s " % argv[0])
        #The command line : ./plotroc.py [-v cv_fold | -T testing_file] [libsvm-options] training_file
        train_file = argv[-1]
        test_file = None
        fold = 5
        options = []
        i = 1
	while i < len(argv)-1:
                if argv[i] == '-T':
                        test_file = argv[i+1]
                        i += 1
                elif argv[i] == '-v':
                        fold = int(argv[i+1])
                        i += 1
                else :
                      	options += [argv[i]]
                i += 1

        return ' '.join(options), fold, train_file, test_file

def plot_roc(deci, label):
        #count of postive and negative labels
        db = []
        pos, neg = 0, 0
        for i in range(len(label)):
                if label[i]>0:
                        pos+=1
                else:
                     	neg+=1
                db.append([deci[i], label[i]])

        #sorting by decision value
        db = sorted(db, key=itemgetter(0), reverse=True)

        #calculate ROC
        xy_arr = []
        tp, fp = 0., 0.                 #assure float division
        for i in range(len(db)):
                if db[i][1]>0:          #positive
                        tp+=1
                else:
                     	fp+=1
                xy_arr.append([fp/neg,tp/pos])

        #area under curve
        aoc = 0.
        prev_x = 0
        for x,y in xy_arr:
                if x != prev_x:
                        aoc += (x - prev_x) * y
                        prev_x = x

        #return aoc, tp/fp pairs (xy_arr)
        return [aoc,xy_arr]



def flattenList(listin,train_file):
        list2=[]
        for item in listin:
            list2.append('\t'.join([str(x) for x in item]+[train_file]))
        final='\n'.join(list2)
        return final


def main():
        if len(argv) <= 1:
                print("Usage: %s [-v cv_fold | -T testing_file] [libsvm-options] training_file" % argv[0])
                raise SystemExit
        param,fold,train_file,test_file = proc_argv()
        output_file = path.split(train_file)[1] + '-rocdata.txt'
        #read data
        train_y, train_x = svm_read_problem(train_file)
        if set(train_y) != set([1,-1]):
                print("ROC is only applicable to binary classes with labels 1, -1")
                raise SystemExit

        #get decision value, with positive = label+
        seed(0) #reset random seed
        if test_file:           #go with test_file
                output_title = "%s on %s" % (path.split(test_file)[1], path.split(train_file)[1])
                test_y, test_x = svm_read_problem(test_file)
                if set(test_y) != set([1,-1]):
                        print("ROC is only applicable to binary classes with labels 1, -1")
                        raise SystemExit
                deci,model = get_pos_deci(train_y, train_x, test_y, test_x, param)
                outdata = plot_roc(deci, test_y)
        else:                           #single file -> CV
                output_title = path.split(train_file)[1]
                deci = get_cv_deci(train_y, train_x, param, fold)
                outdata = plot_roc(deci, train_y)

        print "AUC: %s" % (outdata[0])
        f = open(output_file,'w')
        f.write(flattenList(outdata[1],train_file) +'\n')
        f.close()

if __name__ == '__main__':
        main()
