import numpy as np;

Num_base=3; #number of bases for the oligo
mid=Num_base/2;#The position for the middle base where mutation occurs
N=4**(Num_base-1)*2; #types of mutations

X=np.loadtxt("C:\\MutationSignature\\TCGA\\21bc\\mut.txt")
S=X.shape[0];#number of samples, from input
X=X.reshape((S,N,3))

Genome_BG=np.loadtxt("C:\\MutationSignature\\comp\\base3_simu_sup\\Background.txt")
Genome_BG=Genome_BG.reshape((N,3))[:,0]

p=0.5;#proportion of training data

test_BG=np.random.binomial(Genome_BG.astype(int),1-p)
#Genome_BG-=test_BG;

#temp=test_BG;
#test_BG=Genome_BG;
#Genome_BG=temp;

test_X=np.random.binomial(X.astype(int),1-p)
#X-=test_X;

#temp=test_X;
#test_X=X
#X=temp;

training_error=[];
test_error=[];

for K in range(5,6):
	
	#Hidden variable for EM algorithm
	Z=np.random.random((S,N,3,K));
	Z/=Z.sum(3)[:,:,:,np.newaxis];
	
	#Mutation rate/activity for each signatuer in each sample
	P=np.random.random((S,K))
	
	#Feature matrix
	Feature=np.random.random((Num_base*4-2,K));
	
	def conv_update():
		conv=np.zeros((N,K));
		for i in range(N):
			temp=i;
			for j in range(Num_base):
				if(j<mid):
					conv[i,:]+=Feature[j*4+temp%4,:];
					temp/=4;
				elif(j==mid):
					conv[i,:]+=Feature[j*4+temp%2,:];
					temp/=2;
				elif(j>mid):
					conv[i,:]+=Feature[j*4-2+temp%4,:];
					temp/=4;
		conv=conv+np.abs(conv);
		conv/=2.0;
		return conv;
	
	conv=conv_update();
	
	#signatures
	theta_array=np.random.random((N,K))
	theta_array/=theta_array.sum(0)
	
	#used later for updating
	beta_array=X+10**(-4)
	
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
	
	#Matrix of different types of mutations
	T=np.zeros((N,4*Num_base-2));
	for i in range(N):
		temp=i;
		for j in range(Num_base):
			if(j<mid):
				T[i,4*j+temp%4]=1;
				temp/=4;
			elif(j==mid):
				T[i,4*j+temp%2]=1;
				temp/=2;
			elif(j>mid):
				T[i,4*j-2+temp%4]=1;
				temp/=4;
	
	#conv=np.dot(T,Feature)
	T_inv=np.linalg.pinv(np.dot(np.transpose(T),T));
	Multi_F=np.dot(T_inv,T.transpose());
	Feature=np.dot(Multi_F,theta_array);
	conv=conv_update();			
	
	
	#Regularizer
	for r in range(7):
		
		Reg=10**r-1;
		
		
		
		
		LOSS=((Genome_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
		for i in range(2):
			temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
			temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
			LOSS-=(temp*X[:,type[mid][i],:]).sum()
		
		LOSS+=0.5*Reg*((theta_array-conv)**2).sum()
		
		
		
		
		old=LOSS+1;
		new=LOSS;
		
		
		
		
		while(abs(new-old)>0.01):
			###
			###
			###update signatures
			#one by one for they may take different cycles to converge
			for i in range(K):
				#prepare parameters
				C=Reg*conv[:,i]-Genome_BG*(P[:,i].sum());
				#off-set for thoes equal zero
				temp=np.zeros((S,N,3))
				for j in range(2):
					temp[:,type[mid][j],:]=(P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][j],np.newaxis,:]*Matrix[j,:,:]).sum(3)
				alpha_array=temp
				for j in range(2):
					alpha_array[:,type[mid][j],:]-=P[:,np.newaxis,np.newaxis,i]*theta_array[type[mid][j],np.newaxis,i]*Matrix[j,:,i];
					alpha_array[:,type[mid][j],:]/=P[:,np.newaxis,np.newaxis,i]*Matrix[j,:,i];
				
				#binary search
				#Initialize
				upper=(beta_array/alpha_array).sum((0,2))+C
				v_max=max(upper);
				x_max_max=np.zeros(N);
				x_max_min=np.zeros(N);
				
				v_min=max((beta_array/(alpha_array+1)).sum((0,2))-Reg+C);
				x_min_max=np.ones(N);
				x_min_min=np.zeros(N);
				x_min_min[np.argmax((beta_array/(alpha_array+1)).sum((0,2))-Reg+C)]=1;
				
				#iterate until v_max and v_min are sufficiently close
				while(v_max-v_min>0.0001):
					v_new=1.0/2*(v_max+v_min);
					#v is a convex function of x
					x_new_max=1.0/2*(x_min_max+x_max_max)
					if(not np.isfinite(v_max)):
						v_new=2*np.abs(v_min);
						x_new_max=np.copy(x_min_max)
					x_new_max[v_new>upper]=0 
					x_new_min=np.copy(x_max_min)
					x_new_min[v_new>upper]=0 
					
					#interate x to be determine wether v is smaller or greater than the needed one
					while((x_new_min.sum()-1)*(x_new_max.sum()-1)<0):
						x_new=1.0/2*(x_new_min+x_new_max);
						b_array=(beta_array/(alpha_array+x_new[:,np.newaxis])).sum((0,2))-Reg*x_new+C>v_new
						x_new_min[b_array]=x_new[b_array];
						x_new_max[~b_array]=x_new[~b_array];
					
					if x_new_min.sum()-1>=0 and x_new_max.sum()-1>=0:
						v_min=v_new;
						x_min_min=x_new_min;
						x_min_max=x_new_max;
					else:
						v_max=v_new;
						x_max_min=x_new_min;
						x_max_max=x_new_max;
					
				
				####final update to make x close to the value for v 
				v_new=1.0/2*(v_max+v_min);
				#v is a convex function of x
				x_new_max=1.0/2*(x_min_max+x_max_max)
				x_new_max[v_new>upper]=0
				x_new_min=np.copy(x_max_min)
				x_new_min[v_new>upper]=0
				
				#interate x to be determine wether v is smaller or greater than the needed one
				while(sum((x_new_min-x_new_max)**2)>0.001):
					x_new=1.0/2*(x_new_min+x_new_max);
					b_array=(beta_array/(alpha_array+x_new[:,np.newaxis])).sum((0,2))-Reg*x_new+C>v_new
					x_new_min[b_array]=x_new[b_array];
					x_new_max[~b_array]=x_new[~b_array];
				theta_array[:,i]=1.0/2*(x_new_max+x_new_min)
			
			theta_array/=theta_array.sum(0);
			
			for i in range(2):
				Z[:,type[mid][i],:,:]=P[:,np.newaxis,np.newaxis,:]*theta_array[np.newaxis,type[mid][i],np.newaxis,:]*Matrix[np.newaxis,i,:,:];
			
			Z+=1e-4*Z.sum()/(K*N*3*S);
			Z/=(Z.sum(3))[:,:,:,np.newaxis];
			
			#Update Mutation Matrix
			Matrix=(X[:,type[mid],:,np.newaxis]*Z[:,type[mid],:,:]).sum((0,2));
			Matrix=Matrix+0.01;
			Matrix/=Matrix.sum(1)[:,np.newaxis,:];
			
			#Maximization
			#Update background mutation rate/activity
			P=(X[:,:,:,np.newaxis]*Z).sum((1,2))/(Genome_BG[:,np.newaxis]*theta_array).sum(0);
			
			
			#Update feature
			F_Gradient=np.ones((N,K))
			alpha=1
			while (((alpha*F_Gradient)**2)).sum()>0.01:
				F_Gradient=(theta_array-conv)*np.sign(conv)
				F_Gradient=np.dot(Multi_F,F_Gradient)
				conv_change=np.dot(T,F_Gradient);
				temp=np.dot(T,Feature)/conv_change
				
				alpha=-np.max(temp[temp<0])
				alpha=min(1,alpha*1.0001)
				
				Feature+=alpha*F_Gradient
				conv=conv_update();	
			
			LOSS=((Genome_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
			for i in range(2):
				temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
				temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
				LOSS-=(temp*X[:,type[mid][i],:]).sum()
			
			LOSS+=0.5*Reg*((theta_array-conv)**2).sum()
			
			old=new
			new=LOSS
			
		
		
		print r;
	
	
	
	TEST_LOSS=((test_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
	for i in range(2):
		temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
		temp=np.log(test_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
		TEST_LOSS-=(temp*test_X[:,type[mid][i],:]).sum()
	
	training_error.append(LOSS);
	test_error.append(TEST_LOSS);

np.savetxt("C:\\MutationSignature\\TCGA\\21bc\\ReLU_Feature.txt",Feature, fmt='%s',delimiter='\t')
np.savetxt("C:\\MutationSignature\\TCGA\\21bc\\ReLU_Matrix.txt",Matrix, fmt='%s',delimiter='\t')
np.savetxt("C:\\MutationSignature\\TCGA\\21bc\\ReLU_P.txt",P, fmt='%.18e',delimiter='\t')

#temp=np.copy(Feature[0:3,:])
#Feature[0:3,:]=np.copy(Feature[6:9,:])
#Feature[6:9,:]=temp

conv=conv_update();	
conv_all=np.zeros((N,3,K))
for i in range(2):
	conv_all[type[mid][i]]=conv[type[mid][i],np.newaxis,:]*(Genome_BG[type[mid][i],np.newaxis,np.newaxis]/Genome_BG.sum())*Matrix[i,:,:];

conv_all=conv_all.reshape(N*3,K)
np.savetxt("C:\\MutationSignature\\TCGA\\21bc\\ReLU_conv.txt",conv_all, fmt='%.18e',delimiter='\t')