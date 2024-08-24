# find a SVM line for classification based on the descriptor from SISSO
# read train.dat <and predict.dat>

import numpy as np
import time
from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.model_selection import LeaveOneOut

start=time.time()

# training data
inp=open('train.dat','r').readlines()
X=[]
y_train=np.array([])
for i in range(1,len(inp)):   # from the second line to the last line
	line=str(inp[i]).split()
	feat=[]	
	for j in range(2,len(line)):  # line[0]: name; line[1]: property; line[2] first feature
		feat.append(float(line[j]))
	X.append(feat)
	y_train=np.append(y_train,float(line[1]))
X_train=np.array(X)

# test data
#inp=open('predict.dat','r').readlines()
#X=[]
#y_predict=np.array([])
#for i in range(1,len(inp)):
#        line=str(inp[i]).split()
#        feat=[]
#        for j in range(2,len(line)):
#                feat.append(float(line[j]))
#        X.append(feat)
#        y_predict=np.append(y_predict,float(line[1]))
#X_predict=np.array(X)

# rbf kernel
#clf = GridSearchCV(SVC(kernel='rbf',tol=1e-4,max_iter=-1),
#         cv=10, n_jobs=-1, verbose=1,
#         param_grid={"C": np.logspace(0, 6, 30, base=10), 
#                 "gamma": np.logspace(-6, 2, 30, base=10)},
#         return_train_score=True,refit=True)
# Linear kernel
clf = GridSearchCV(SVC(kernel='linear',tol=1e-4,max_iter=-1),
         cv=10, n_jobs=-1, verbose=1,
         param_grid={"C": np.logspace(1, 5, 10, base=10)},
         return_train_score=True,refit=True)

print("training:")
clf.fit(X_train, y_train)  
print("params: ")
print(clf.cv_results_['params'])  
print("mean_test_score: ")
print(clf.cv_results_['mean_test_score']) 
print("mean_train_score: ")
print(clf.cv_results_['mean_train_score']) 
print("best CV score: ")
print(clf.best_score_) 
print("best_params: ")
print(clf.best_params_) 
print("best_estimator: ")
print(clf.best_estimator_) 
print("coef: (sum(c_i*x_i)+c_0=0) ")
print(clf.best_estimator_.coef_)
print("intercept: ")
print(clf.best_estimator_.intercept_)
print("support: ")
print(clf.best_estimator_.support_)
print(clf.best_estimator_.support_vectors_)

def train_output(X_train,y_train):
    y_pred=clf.predict(X_train)
    er=0
    fail_index=[]
    for i in range(len(y_train)):
          print('Y_ture, Y_pred, data_point: ',np.around(y_train[i],decimals=3),np.around(y_pred[i],decimals=3),i)
          if np.sign(y_train[i]) != np.sign(y_pred[i]):
              er+=1
              fail_index.append(i)
    print('Wrong prediction: Data point ',fail_index)
    print("Total number of misclassified data: ",er)
train_output(X_train,y_train)
print("training_score (accuracy): ",clf.best_estimator_.score(X_train,y_train))


#def validate_output(X_predict,y_predict):
#    y_pred=clf.predict(X_predict)
#    er=0
#    fail_index=[]
#    for i in range(len(y_predict)):
#          print('Y_ture, Y_pred, data_point: ',np.around(y_predict[i],decimals=3),np.around(y_pred[i],decimals=3),i)
#          if np.sign(y_predict[i]) != np.sign(y_pred[i]):
#              er+=1
#              fail_index.append(i)
#    print('Wrong prediction: Data point ',fail_index)
#    print("Total number of misclassified data: ",er)

#print("validation: ")
#validate_output(X_predict,y_predict)
#print("validation_score (accuracy): ",clf.best_estimator_.score(X_predict,y_predict))



end=time.time()
print("wall-clock time (seconds):",end-start)
#
