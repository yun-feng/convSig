
import numpy as np;

Num_base=3; #number of bases for the oligo
mid=Num_base/2;#The position for the middle base where mutation occurs
N=4**(Num_base-1)*2; #types of mutations

X=np.loadtxt("C:\\MutationSignature\\TCGA\\21bc\\mut.txt")
S=X.shape[0];#number of samples, from input
X=X.reshape((S,N,3))


Genome_BG=np.loadtxt("C:\\MutationSignature\\comp\\base3_simu_sup_CRC\\Background.txt")
#Genome_BG=Genome_BG.reshape((N,3))[:,0]

p=0.5;#proportion of training data

test_BG=np.random.binomial(Genome_BG.astype(int),1-p)
#Genome_BG-=test_BG;

#temp=test_BG;
#test_BG=Genome_BG;
#Genome_BG=temp;

test_X=np.random.binomial(X.astype(int),1-p)
#X-=test_X;

X+=10^(-5)*X.sum()/np.size(X)
test_X+=10^(-5)*test_X.sum()/np.size(X)

#temp=test_X;
#test_X=X
#X=temp;

training_error=[];
test_error=[];

for K in range(5,6):
	
	#K=6;#number of filters, decided by cross-validation
	
	
	#Hidden variable for EM algorithm
	Z=np.random.random((S,N,3,K));
	Z/=Z.sum(3)[:,:,:,np.newaxis];
	
	#Mutation rate/activity for each signatuer in each sample
	P=np.random.random((S,K))
	
	#Laugrange Multiplier
	#Mu=np.random.random(S);
	
	#Define Mutation Siganture/Filter
	
	#Feature matrix
	Feature=np.random.random((Num_base,4,K));
	#The middle base only has two probabilitys
	Feature[mid,2:4,:]=0;
	#entry is converted to the exponential form:eta
	#each base position sums to 1
	Feature=Feature/Feature.sum(1)[:,np.newaxis,:]
	
	#pre-compute the convolution results after transformation
	def conv_update():
		conv=np.ones((N,K));
		for i in range(N):
			temp=i;
			for j in range(Num_base):
				if(j==mid):
					conv[i,:]*=Feature[j,temp%2,:];
					temp/=2;
				else:
					conv[i,:]*=Feature[j,temp%4,:];
					temp/=4;
		return conv;
	
	conv=conv_update();
	
	#Mutation matrix
	#each base position sums to 1
	Matrix=np.random.random((2,3,K))
	Matrix/=Matrix.sum(1)[:,np.newaxis,:]
	
	
	#store the index for each type of mutation to each particular pattern
	type=[];
	base_div=1;
	for i in range(Num_base):
		if(i==mid):
			temp=[[],[]];
			base_mod=2;
		else:
			temp=[[],[],[],[]];
			base_mod=4;
			
		for j in range(N):
			base_type=j/(base_div) % base_mod;
			temp[base_type].append(j);
		
		if(i==mid):
			base_div*=2;
		else:
			base_div*=4;
			
		type.append(temp);	
	
	
	LOSS=((Genome_BG[:,np.newaxis]*conv).sum(0)*P).sum()
	for i in range(2):
		temp=P[:,np.newaxis,np.newaxis,:]*conv[np.newaxis,type[mid][i],np.newaxis,:]*Matrix[np.newaxis,i,:,:]
		temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
		LOSS-=(temp*X[:,type[mid][i],:]).sum()
	
	complete_loss=((Genome_BG[:,np.newaxis]*conv).sum(0)*P).sum()
	complete_loss+=((Z*np.log(Z)).sum(3)*X).sum()
	for i in range(2):
		temp=P[:,np.newaxis,np.newaxis,:]*conv[np.newaxis,type[mid][i],np.newaxis,:]*Matrix[np.newaxis,i,:,:]
		complete_loss-=((Z[:,type[mid][i],:,:]*np.log(temp)).sum(3)*X[:,type[mid][i],:]).sum()
	
	old=0;
	new=LOSS;
	LOSS_array=[];	
	LOSS_array.append(LOSS);
	
	#test error
	test_LOSS=((test_BG[:,np.newaxis]*conv).sum(0)*P).sum()
	for i in range(2):
		temp=P[:,np.newaxis,np.newaxis,:]*conv[np.newaxis,type[mid][i],np.newaxis,:]*Matrix[np.newaxis,i,:,:]
		temp=np.log(test_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
		test_LOSS-=(temp*test_X[:,type[mid][i],:]).sum()	
	
	test_array=[test_LOSS];
	#	
	#iterate until convergence
	#
	while(abs(new-old)>10**(-3)):
		#Expectation Step
		for i in range(2):
			Z[:,type[mid][i],:,:]=P[:,np.newaxis,np.newaxis,:]*conv[np.newaxis,type[mid][i],np.newaxis,:]*Matrix[np.newaxis,i,:,:];
		
		Z/=(Z.sum(3))[:,:,:,np.newaxis];
		
		
		
		#Update Mutation Matrix
		Matrix=(X[:,type[mid],:,np.newaxis]*Z[:,type[mid],:,:]).sum((0,2))
		
		Matrix/=Matrix.sum(1)[:,np.newaxis,:]
		
		#Update Feature Matrix
		Product=(X[:,:,:,np.newaxis]*Z[:,:,:,:]).sum((0,2))
		for i in range(Num_base):
			if(i==mid):
				Feature[i,0:2,:]=1
			else:
				Feature[i,:,:]=1
			conv=conv_update();
			
			Weight=P.sum(0)[np.newaxis,:]*Genome_BG[:,np.newaxis]*conv	
			Weight_i=Weight[type[i],:].sum(1)
				
			Product_i=Product[type[i],:].sum(1)
			
			if(i==mid):
				Feature[i,0:2,:]=Product_i/Weight_i;
			else:
				Feature[i,:,:]=Product_i/Weight_i;
			Feature[i,:,:]/=Feature[i].sum(0)[np.newaxis,:];
			
			
			
		conv=conv_update();
		
		#Maximization
		#Update background mutation rate/activity
		P=(X[:,:,:,np.newaxis]*Z).sum((1,2))/(Genome_BG[:,np.newaxis]*conv).sum(0)
		
		
		old=LOSS
		#Loss Function Evaluation
		LOSS=((Genome_BG[:,np.newaxis]*conv).sum(0)*P).sum()
		for i in range(2):
			temp=P[:,np.newaxis,np.newaxis,:]*conv[np.newaxis,type[mid][i],np.newaxis,:]*Matrix[np.newaxis,i,:,:]
			temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
			LOSS-=(temp*X[:,type[mid][i],:]).sum()
		
		new=LOSS;
		if(new>old):
			print "Something Wrong!";
		LOSS_array.append(LOSS);
		
		#test error
		test_LOSS=((test_BG[:,np.newaxis]*conv).sum(0)*P).sum()
		for i in range(2):
			temp=P[:,np.newaxis,np.newaxis,:]*conv[np.newaxis,type[mid][i],np.newaxis,:]*Matrix[np.newaxis,i,:,:]
			temp=np.log(test_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
			test_LOSS-=(temp*test_X[:,type[mid][i],:]).sum()
		test_array.append(test_LOSS);
	
	training_error.append(LOSS);
	test_error.append(test_LOSS);
	print K;


np.savetxt("C:\\MutationSignature\\TCGA\\21bc\\exp_Feature.txt",Feature, fmt='%s',delimiter='\t')
np.savetxt("C:\\MutationSignature\\TCGA\\21bc\\exp_Matrix.txt",Matrix, fmt='%s',delimiter='\t')
np.savetxt("C:\\MutationSignature\\TCGA\\21bc\\exp_p.txt",P, fmt='%.18e',delimiter='\t')

def conv_reverse():
	conv=np.ones((N,K));
	for i in range(N):
		temp=i;
		for j in range(Num_base):
			if(j==mid):
				conv[i,:]*=Feature[j,temp%2,:];
				temp/=2;
			else:
				conv[i,:]*=Feature[Num_base-j-1,temp%4,:];
				temp/=4;
	return conv;

#conv=conv_reverse();
conv_all=np.zeros((N,3,K))
for i in range(2):
	conv_all[type[mid][i]]=conv[type[mid][i],np.newaxis,:]*(Genome_BG[type[mid][i],np.newaxis,np.newaxis]/Genome_BG.sum())*Matrix[i,:,:];

conv_all=conv_all.reshape(N*3,K)
np.savetxt("C:\\MutationSignature\\TCGA\\21bc\\exp_conv.txt",conv_all, fmt='%.18e',delimiter='\t')
