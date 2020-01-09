import argparse
import random
import scipy.interpolate
import numpy as np
from sklearn.linear_model import LogisticRegression

def minimize_error(positives,negatives):
	if positives == 0: return False,"NA","NA",0 # a region with no positives has zero classification power
	if negatives == 0: negatives = [random.choice(range(0,50)) for i in range(len(positives))]
	pos_min = min(positives)
	neg_max = max(negatives)
	dt = float(pos_min+neg_max)/2.0
	positives = [i-dt for i in positives]
	negatives = [i-dt for i in negatives]
	x_train = np.array(positives+negatives).reshape(-1,1)
	y_train = np.array([1 for i in range(len(positives))]+[0 for i in range(len(negatives))])
	LR = LogisticRegression(class_weight='balanced').fit(x_train,y_train)
	classification_power = LR.score(x_train,y_train)
	slope = float(LR.coef_[0][0])
	y_int = float(LR.intercept_[0])
	return(dt,slope,y_int,classification_power)


def process_results(results):
	if results == {}: return None
	else:
		with open(args.o,'a') as out:
			for uniref,x in sorted(results.items()):
				spline = []
				for start,y in x.items():
					neg,pos = 0,0
					for label,slopes in y.items():
						if 'N' in label:
							neg = slopes
						else:
							pos = slopes
					decision_thresh,lr_slope,lr_y_int,lr_classification_power = minimize_error(pos,neg)
					if decision_thresh != False:
						spline.append([start,decision_thresh,lr_slope,lr_y_int,lr_classification_power])
				if spline != []:
					spline = sorted(spline)
					spline_starts,spline_dt,LR_slopes,LR_y_ints,LR_classification_power = zip(*spline)
					if len(spline_starts) > 1:
						spline_model = scipy.interpolate.CubicSpline(spline_starts, spline_dt, bc_type='natural')
						# cubic polynomial is represented as ax^3 + bx^2 + cx + d
						a = spline_model.c[0]
						b = spline_model.c[1]
						c = spline_model.c[2]
						d = spline_model.c[3]
						out.write(uniref+'\t'+','.join(map(str,spline_model.x))+'\t'+','.join(map(str,a))+'\t'+','.join(map(str,b))+'\t'+','.join(map(str,c))+'\t'+','.join(map(str,d))+'\t'+','.join(map(str,LR_slopes))+'\t'+','.join(map(str,LR_y_ints))+'\t'+','.join(map(str,LR_classification_power))+'\n')


parser = argparse.ArgumentParser()
parser.add_argument("-s",help="input slope results file from training data")
parser.add_argument("-o",help="desired output file name")
args = parser.parse_args()

with open(args.o,'a') as out:
	out.write("Marker_gene_group\tKnots\tA\tB\tC\tD\tLR_slopes\tLR_y_ints\tClassification_power\n")

results ={}

# Example of args.s file structure
# EOG092C00SU_ACAROMYCES  1060    1.046875	NG
for i in open(args.s):
	tmp = i.strip().split('\t')
	if len(tmp) < 4: continue
	uniref = tmp[0]
	if uniref == '0': continue
	start = int(tmp[1]) - int(tmp[1])%10
	slope = float(tmp[2])
	label = tmp[3]
	if 'N' in label: label = 'N'
	if uniref not in results:
		process_results(results)
		results = {}
		results[uniref] = {start:{label:[slope]}}
	else:
		if start not in results[uniref]:
			results[uniref][start] = {label:[slope]}
		else:
			if label not in results[uniref][start]:
				results[uniref][start][label] = [slope]
			else:
				results[uniref][start][label].append(slope)

process_results(results)



