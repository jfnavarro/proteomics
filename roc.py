from pylab import *
import random

def _remove_duplicate_styles(rocList):
    
    pref_styles = ['cx-','mx-','yx-','gx-','bx-','rx-']
    points = 'ov^>+xd'
    colors = 'bgrcmy'
    lines = ['-','-.',':']
    
    rand_ls = []
    
    for r in rocList:
        if r.linestyle not in rand_ls:
            rand_ls.append(r.linestyle)
        else:
            while True:
                if len(pref_styles) > 0:
                    pstyle = pref_styles.pop()
                    if pstyle not in rand_ls:
                        r.linestyle = pstyle
                        rand_ls.append(pstyle)
                        break
                else:
                    ls = ''.join(random.sample(colors,1) + random.sample(points,1)+ random.sample(lines,1))
                    if ls not in rand_ls:
                        r.linestyle = ls
                        rand_ls.append(ls)
                        break

def plot_multiple_roc(rocList,titlep='',labels=None, include_baseline=False, equal_aspect=True,filename=""):
    """ Plots multiple ROC curves on the same chart. 
        Parameters:
            rocList: the list of ROCData objects
            title: The tile of the chart
            labels: The labels of each ROC curve
            include_baseline: if it's  True include the random baseline
            equal_aspect: keep equal aspect for all roc curves
    """
    clf()
    ylim((0,1))
    xlim((0,1))
    xticks(arange(0,1.1,.1))
    yticks(arange(0,1.1,.1))
    grid(True)
    if equal_aspect:
        cax = gca()
        cax.set_aspect('equal')
    xlabel("1 - Specificity")
    ylabel("Sensitivity")
    title(titlep)
    if not labels:
        labels = [ '' for x in rocList]
    _remove_duplicate_styles(rocList)
    for ix, r in enumerate(rocList):
        plot([x[0] for x in r.derived_points], [y[1] for y in r.derived_points], r.linestyle, linewidth=1, label=labels[ix])
    if include_baseline:
        plot([0.0,1.0], [0.0, 1.0], 'k-', label= 'random')
    if labels:
        legend(loc='lower right')
    
    if(filename != ""):
        savefig(filename, format='png')
    else:
        show()
    
    
    
class ROCData(object):
    """ Class that generates an ROC Curve for the data.
        Data is in the following format: a list l of tutples t
        where:
            t[0] = 1 for positive class and t[0] = 0 for negative class
            t[1] = score
            t[2] = label
    """
    def __init__(self,data,linestyle='rx-'):
        """ Constructor takes the data and the line style for plotting the ROC Curve.
            Parameters:
                data: The data a listl of tuples t (l = [t_0,t_1,...t_n]) where:
                      t[0] = 1 for positive class and 0 for negative class
                      t[1] = a score
                       t[2] = any label (optional)
                lineStyle: THe matplotlib style string for plots.
                
            Note: The ROCData is still usable w/o matplotlib. The AUC is still available, 
                  but plots cannot be generated.
        """
        self.data = sorted(data,lambda x,y: cmp(y[1],x[1]))
        self.linestyle = linestyle
        self.auc() #Seed initial points with default full ROC
    
    def auc(self,fpnum=0):
        """ Uses the trapezoidal ruel to calculate the area under the curve. If fpnum is supplied, it will 
            calculate a partial AUC, up to the number of false positives in fpnum (the partial AUC is scaled
            to between 0 and 1).
            It assumes that the positive class is expected to have the higher of the scores (s(+) < s(-))
            Parameters:
                fpnum: The cumulativr FP count (fps)
            Return:
            
        """
        fps_count = 0
        relevant_pauc = []
        current_index = 0
        max_n = len([x for x in self.data if x[0] == 0])
        if fpnum == 0:
            relevant_pauc = [x for x in self.data]
        elif fpnum > max_n:
            fpnum = max_n
        #Find the upper limit of the data that does not exceed n FPs
        else:
            while fps_count < fpnum:
                relevant_pauc.append(self.data[current_index])
                if self.data[current_index][0] == 0:
                    fps_count += 1
                current_index +=1
        total_n = len([x for x in relevant_pauc if x[0] == 0])
        total_p = len(relevant_pauc) - total_n
        
        #Convert to points in a ROC
        previous_df = -1000000.0
        current_index = 0
        points = []
        tp_count, fp_count = 0.0 , 0.0
        tpr, fpr = 0, 0
        while current_index < len(relevant_pauc):
            df = relevant_pauc[current_index][1]
            if previous_df != df:
                points.append((fpr,tpr,fp_count))
            if relevant_pauc[current_index][0] == 0:
                fp_count +=1
            elif relevant_pauc[current_index][0] == 1:
                tp_count +=1
            fpr = fp_count/total_n
            tpr = tp_count/total_p
            previous_df = df
            current_index +=1
        points.append((fpr,tpr,fp_count)) #Add last point
        points.sort(key=lambda i: (i[0],i[1]))
        self.derived_points = points
        
        return self._trapezoidal_rule(points)


    def _trapezoidal_rule(self,curve_pts):
        """ Method to calculate the area under the ROC curve"""
        cum_area = 0.0
        for ix,x in enumerate(curve_pts[0:-1]):
            cur_pt = x
            next_pt = curve_pts[ix+1]
            cum_area += ((cur_pt[1]+next_pt[1])/2.0) * (next_pt[0]-cur_pt[0])
        return cum_area
        
    def calculateStandardError(self,fpnum=0):
        """ Returns the standard error associated with the curve.
            Parameters:
                fpnum: The cumulativr FP count (fps)
            Return:
                the standard error.
        """
        area = self.auc(fpnum)
        
        #real positive cases
        Na =  len([ x for x in self.data if x[0] == 1])
        
        #real negative cases
        Nn =  len([ x for x in self.data if x[0] == 0])
        
        
        Q1 = area / (2.0 - area)
        Q2 = 2 * area * area / (1.0 + area)
        
        return math.sqrt( ( area * (1.0 - area)  +   (Na - 1.0) * (Q1 - area*area) +
                        (Nn - 1.0) * (Q2 - area * area)) / (Na * Nn))
                            
    
    def plot(self,titlep='',include_baseline=False,equal_aspect=True):
        """ Method that generates a plot of the ROC curve 
            Parameters:
                title: Title of the chart
                include_baseline: Add the baseline plot line if it's True
                equal_aspect: Aspects to be equal for all plot
        """
        
        clf()
        plot([x[0] for x in self.derived_points], [y[1] for y in self.derived_points], self.linestyle)
        if include_baseline:
            plot([0.0,1.0], [0.0,1.0],'k-.')
        ylim((0,1))
        xlim((0,1))
        xticks(arange(0,1.1,.1))
        yticks(arange(0,1.1,.1))
        grid(True)
        if equal_aspect:
            cax = gca()
            cax.set_aspect('equal')
        xlabel('1 - Specificity')
        ylabel('Sensitivity')
        title(titlep)
        
        show()
        
    
    def confusion_matrix(self,threshold,do_print=False):
        """ Returns the confusion matrix (in dictionary form) for a fiven threshold
            where all elements > threshold are considered 1 , all else 0.
            Parameters:
                threshold: threshold to check the decision function
                do_print:  if it's True show the confusion matrix in the screen
            Return:
                the dictionary with the TP, FP, FN, TN
        """
        pos_points = [x for x in self.data if x[1] >= threshold]
        neg_points = [x for x in self.data if x[1] < threshold]
        tp,fp,fn,tn = self._calculate_counts(pos_points,neg_points)
        if do_print:
            print "\t Actual class"
            print "\t+(1)\t-(0)"
            print "+(1)\t%i\t%i\tPredicted" % (tp,fp)
            print "-(0)\t%i\t%i\tclass" % (fn,tn)
        return {'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn}
        

    
    def evaluateMetrics(self,matrix,metric=None,do_print=False):
        """ Returns the metrics evaluated from the confusion matrix.
            Parameters:
                matrix: the confusion matrix
                metric: the specific metric of the default value is None (all metrics).
                do_print:  if it's True show the metrics in the screen
            Return:
                the dictionary with the Accuracy, Sensitivity, Specificity,Efficiency,
                                        PositivePredictiveValue, NegativePredictiveValue, PhiCoefficient
        """
        
        accuracy = (matrix['TP'] + matrix['TN'])/ float(sum(matrix.values()))
        
        sensitivity = (matrix['TP'])/ float(matrix['TP'] + matrix['FN'])
        
        specificity = (matrix['TN'])/float(matrix['TN'] + matrix['FP'])
        
        efficiency = (sensitivity + specificity) / 2.0
        
        positivePredictiveValue =  matrix['TP'] / float(matrix['TP'] + matrix['FP'])

        NegativePredictiveValue = matrix['TN'] / float(matrix['TN'] + matrix['FN'])
        
        PhiCoefficient = (matrix['TP'] * matrix['TN'] - matrix['FP'] * matrix['FN'])/(
                            math.sqrt( (matrix['TP'] + matrix['FP']) *
                                       (matrix['TP'] + matrix['FN']) *
                                       (matrix['TN'] + matrix['FP']) *
                                       (matrix['TN'] + matrix['FN']))) or 1.0
                                    
        if do_print:
            print 'Sensitivity: ' , sensitivity
            print 'Specificity: ' , specificity
            print 'Efficiency: ' , efficiency
            print 'Accuracy: ' , accuracy
            print 'PositivePredictiveValue: ' , positivePredictiveValue
            print 'NegativePredictiveValue' , NegativePredictiveValue
            print 'PhiCoefficient' , PhiCoefficient
            
        
        return {'SENS': sensitivity, 'SPEC': specificity, 'ACC': accuracy, 'EFF': efficiency,
                'PPV':positivePredictiveValue, 'NPV':NegativePredictiveValue , 'PHI':  PhiCoefficient}


    def _calculate_counts(self,pos_data,neg_data):
        """ Calculates the number of false positives, true positives, false negatives and true negatives """
        tp_count = len([x for x in pos_data if x[0] == 1])
        fp_count = len([x for x in pos_data if x[0] == 0])
        fn_count = len([x for x in neg_data if x[0] == 1])
        tn_count = len([x for x in neg_data if x[0] == 0])
        return tp_count,fp_count,fn_count, tn_count
        