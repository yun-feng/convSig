import numpy as np;

Num_base=5; #number of bases for the oligo
mid=Num_base/2;#The position for the middle base where mutation occurs
N=4**Num_base; #types of mutations

X=np.loadtxt("C:\\MutationSignature\\comp\\base5_conv\\Simu3.txt")

S=X.shape[0];#number of samples, from input
X=X.reshape((S,N,3))

Genome_BG=np.loadtxt("C:\\MutationSignature\\comp\\base5_conv\\Background.txt")
Genome_BG=Genome_BG.reshape((N,3))[:,0]


##special matrice
temp=np.array([0]*Num_base*4*(Num_base-2)*4);
temp[map(lambda x:x*(Num_base-2)*4+x,range((Num_base-2)*4))]=1
I1=temp;
I1.resize(Num_base*4,(Num_base-2)*4);

temp=np.array([0]*Num_base*4*(Num_base-2)*4);
temp[map(lambda x:(x+4)*(Num_base-2)*4+x,range((Num_base-2)*4))]=1
I2=temp;
I2.resize(Num_base*4,(Num_base-2)*4);

temp=np.array([0]*Num_base*4*(Num_base-2)*4);
temp[map(lambda x:(x+8)*(Num_base-2)*4+x,range((Num_base-2)*4))]=1
I3=temp;
I3.resize(Num_base*4,(Num_base-2)*4);


T=np.zeros((N,4*Num_base));

for i in range(N):
	temp=i;
	for j in range(Num_base):
		T[i,4*j+temp%4]=1;
		temp/=4;

#store the index for each type of mutation to each particular pattern
type=[];
base_div=1;
for i in range(Num_base):
	temp=[[],[],[],[]];
	
	for j in range(N):
		base_type=j/(base_div) % 4;
		temp[base_type].append(j);
	
	base_div*=4;
	
	type.append(temp);	

#Restict Ratio(0.5~1)

Ratio=0.6;


for K in range(5,6):
	
	#Hidden variable for EM algorithm
	Z=np.random.random((S,N,3,K));
	Z/=Z.sum(3)[:,:,:,np.newaxis];
	
	#Mutation rate/activity for each signatuer in each sample
	P=np.random.random((S,K))
	
	#Feature matrix
	Feature=np.ones(((Num_base-2)*4,K));
	Layer2=np.ones((2,K));
	A_list=Layer2[0]*I1[:,:,np.newaxis]+I2[:,:,np.newaxis]+Layer2[1]*I3[:,:,np.newaxis];
	
	conv=np.zeros((N,K));
	
	for i in range(K):
		conv[:,i]=np.dot(T,np.dot(A_list[:,:,i],Feature[:,i]));
	
	conv=np.exp(conv);
	
	#signatures
	theta_array=np.random.random((N,K))
	theta_array/=theta_array.sum(0)
	
	#used later for updating
	beta_array=X+10**(-4)
	
	#Mutation matrix
	#each base position sums to 1
	Matrix=np.random.random((4,3,K))
	Matrix/=Matrix.sum(1)[:,np.newaxis,:]
	
	#Regularizer
	for r in range(7):
		
		Reg=10**r-1;
		
		
		LOSS=((Genome_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
		for i in range(4):
			temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
			temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
			LOSS-=(temp*X[:,type[mid][i],:]).sum()
		
		LOSS+=0.5*Reg*((theta_array-conv)**2).sum()
		
		old=LOSS+2;
		new=LOSS;
		
		
		
		
		while(abs(new-old)>0.01):
			###
			###
			###update signatures
			#one by one for they may take different cycles to converge
			old=new
			if np.any(conv<theta_array*Ratio):
				print "wrong";
			
			
			
			for i in range(K):
				#prepare parameters
				C=Reg*conv[:,i]-Genome_BG*(P[:,i].sum());
				#off-set for thoes equal zero
				temp=np.zeros((S,N,3))
				for j in range(4):
					temp[:,type[mid][j],:]=(P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][j],np.newaxis,:]*Matrix[j,:,:]).sum(3)
				alpha_array=temp
				for j in range(4):
					alpha_array[:,type[mid][j],:]-=P[:,np.newaxis,np.newaxis,i]*theta_array[type[mid][j],np.newaxis,i]*Matrix[j,:,i];
					alpha_array[:,type[mid][j],:]/=P[:,np.newaxis,np.newaxis,i]*Matrix[j,:,i];
				
				#binary search
				#Initialize
				upper=(beta_array/alpha_array).sum((0,2))+C
				lower=(beta_array/(alpha_array+1/Ratio*(conv[:,i])[:,np.newaxis])).sum((0,2))-Reg*1/Ratio*(conv[:,i])+C
				
				v_max=max(upper);
				x_max_max=np.zeros(N);
				x_max_min=np.zeros(N);
				
				v_min=min(lower)
				x_min_max=1.0/Ratio*conv[:,i];
				x_min_min=1.0/Ratio*conv[:,i];
				
				
				#iterate until v_max and v_min are sufficiently close
				while(v_max-v_min>=1e-6*(np.abs(v_max)+np.abs(v_min))):
					v_new=1.0/2*(v_max+v_min);
					#v is a convex function of x
					x_new_max=1.0/2*(x_min_max+x_max_max)
					if(not np.isfinite(v_max)):
						v_new=2*np.abs(v_min);
						x_new_max=np.copy(x_min_max)
					
					x_new_min=np.copy(x_max_min)
					if(np.any(v_new>upper)):
						x_new_max[v_new>upper]=0 
						x_new_min[v_new>upper]=0 
					if (np.any(v_new<lower)):
						x_new_max[v_new<lower]=1.0/Ratio*conv[v_new<lower,i]
						x_new_min[v_new<lower]=1.0/Ratio*conv[v_new<lower,i]
					
					
					
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
				x_new_min=np.copy(x_max_min)
				
				#interate x to be determine wether v is smaller or greater than the needed one
				while(sum((x_new_min-x_new_max)**2)>1e-6*sum((np.abs(x_new_min)+np.abs(x_new_max))**2)):
					x_new=1.0/2*(x_new_min+x_new_max);
					b_array=(beta_array/(alpha_array+x_new[:,np.newaxis])).sum((0,2))-Reg*x_new+C>v_new
					x_new_min[b_array]=x_new[b_array];
					x_new_max[~b_array]=x_new[~b_array];
				theta_array[:,i]=1.0/2*(x_new_max+x_new_min)
				theta_array[:,i]/=theta_array[:,i].sum()
			
			#LOSS=((Genome_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
			#for i in range(4):
			#	temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
			#	temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
			#	LOSS-=(temp*X[:,type[mid][i],:]).sum()
			
			#LOSS+=0.5*Reg*((theta_array-conv)**2).sum()
			
			#old1=new
			#new=LOSS
			#if new>old1:
			#	print "wrong",new,old1;
			
			#impose nonzero	
			
			if (np.any(np.log(conv)-np.log(Ratio*theta_array)<=0)):
				theta_array[np.log(conv)-np.log(Ratio*theta_array)<=0]=1.0/((1.0+1e-5)*Ratio)*conv[np.log(conv)-np.log(Ratio*theta_array)<=0]
			
			for i in range(4):
				Z[:,type[mid][i],:,:]=P[:,np.newaxis,np.newaxis,:]*theta_array[np.newaxis,type[mid][i],np.newaxis,:]*Matrix[np.newaxis,i,:,:];
			Z+=1e-4*Z.sum()/(K*N*3*S);
			Z/=(Z.sum(3))[:,:,:,np.newaxis];
			
			#print "theta";
			
			
			
			
			
			#Update Mutation Matrix
			Matrix=(X[:,type[mid],:,np.newaxis]*Z[:,type[mid],:,:]).sum((0,2));
			Matrix=Matrix+0.01;
			Matrix/=Matrix.sum(1)[:,np.newaxis,:];
			
			#Maximization
			#Update background mutation rate/activity
			P=(X[:,:,:,np.newaxis]*Z).sum((1,2))/(Genome_BG[:,np.newaxis]*theta_array).sum(0);
			
			#LOSS=((Genome_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
			#for i in range(4):
			#	temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
			#	temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
			#	LOSS-=(temp*X[:,type[mid][i],:]).sum()
			
			#LOSS+=0.5*Reg*((theta_array-conv)**2).sum()
			
			#old1=new
			#new=LOSS
			#if new>old1:
			#	print "wrong",new,old1;
				
			
			for i in range(K):
				#Update feature
				A=A_list[:,:,i];
				
				AT=np.dot(T,A).transpose()
				
				#ignore where theta equals zero
				Include_array=theta_array[:,i]>0
				
				for Neu in range(18):
					
					Neu_mul=1.0/(10**Neu);
					
					N_loss=0.5*((theta_array[:,i]-conv[:,i])**2).sum()-Neu_mul*np.log(np.dot(Feature[:,i],AT)[Include_array]-np.log(Ratio*theta_array[Include_array,i])).sum()
					
					F_Gradient=np.ones(N)
					
					while True:
						
						F_Gradient=(-(theta_array[:,i]-conv[:,i])*conv[:,i]-Neu_mul/(np.dot(Feature[:,i],AT)-np.log(Ratio*theta_array[:,i])))*AT
						
						Neu_Div=2*conv[:,i]**2-theta_array[:,i]*conv[:,i]+Neu_mul/((np.dot(Feature[:,i],AT)-np.log(Ratio*theta_array[:,i]))**2)
						Hessian=np.dot(AT*Neu_Div,AT.transpose())
						Hessian=np.linalg.pinv(Hessian);
						
						F_Gradient=F_Gradient.sum(1)
						
						Neu_F=np.dot(Hessian,F_Gradient)
						N_dec=np.dot(F_Gradient,Neu_F)
						
						if (Neu_F**2).sum()/(Feature[:,i]**2).sum()<1e-4:
							break;
						
						alpha=1;
						
						while True:
							Temp_F=Feature[:,i]-alpha*Neu_F;
							conv[:,i]=np.dot(T,np.dot(A,Temp_F));
							conv[:,i]=np.exp(conv[:,i]);
							
							New_N_loss=0.5*((theta_array[:,i]-conv[:,i])**2).sum()-Neu_mul*np.log(np.dot(Temp_F,AT)[Include_array]-np.log(Ratio*theta_array[Include_array,i])).sum();
							if N_loss-New_N_loss>=0.4*alpha*N_dec:
								break;
							
							alpha=alpha*0.8;
						
						N_loss=New_N_loss;
						Feature[:,i]=Temp_F;
					
					
						
						
			#print "Feature";	
			#LOSS=((Genome_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
			#for i in range(4):
			#	temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
			#	temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
			#	LOSS-=(temp*X[:,type[mid][i],:]).sum()
			
			#LOSS+=0.5*Reg*((theta_array-conv)**2).sum()
			
			#old1=new
			#new=LOSS
			#if new>old1:
			#	print "wrong",new,old1;			
			##update Layer2
			##I1
			
			
			for i in range(K):
				a_array=np.dot(T,np.dot(I1,Feature[:,i]))
				b_array=np.dot(T,np.dot(I2,Feature[:,i]))+Layer2[1,i]*np.dot(T,np.dot(I3,Feature[:,i]))
				temp=-b_array+np.log(Ratio*theta_array[:,i]);
				temp=temp/a_array;
				x_max=np.inf;
				x_min=-np.inf
				if(np.any(a_array>0)):
					x_min=np.max(temp[a_array>0]);
				if(np.any(a_array<0)):
					x_max=np.min(temp[a_array<0]);
				
				while not x_max-x_min<1e-7*(np.abs(x_max)+np.abs(x_min)):
					x_new=0.5*(x_max+x_min);
					if (not np.isfinite(x_max)):
						x_new=2*np.abs(x_min+1);
					if (not np.isfinite(x_min)):
						x_new=-2*np.abs(x_max+1);
					Grad_x=(np.exp(x_new*a_array+b_array)-theta_array[:,i])*np.exp(x_new*a_array+b_array)*a_array;
					Grad_x=Grad_x.sum()
					if(Grad_x<0):
						x_min=x_new;
					else:
						x_max=x_new;
				
				Layer2[0,i]=0.5*(x_max+x_min);
				#A_list[:,:,i]=Layer2[0,i]*I1+I2+Layer2[1,i]*I3
				
				#conv[:,i]=np.dot(T,np.dot(A_list[:,:,i],Feature[:,i]));
	
				#conv[:,i]=np.exp(conv[:,i]);
			
			#print "Layer2_0";
			
			#LOSS=((Genome_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
			#for i in range(4):
			#	temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
			#	temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
			#	LOSS-=(temp*X[:,type[mid][i],:]).sum()
			
			#LOSS+=0.5*Reg*((theta_array-conv)**2).sum()
			
			#old1=new
			#new=LOSS
			#if new>old1:
			#	print "wrong",new,old1;
			##I3
			
			
			for i in range(K):
				a_array=np.dot(T,np.dot(I3,Feature[:,i]))
				b_array=np.dot(T,np.dot(I2,Feature[:,i]))+Layer2[0,i]*np.dot(T,np.dot(I1,Feature[:,i]))
				temp=-b_array+np.log(Ratio*theta_array[:,i]);
				temp=temp/a_array;
				x_max=np.inf;
				x_min=-np.inf
				if(np.any(a_array>0)):
					x_min=np.max(temp[a_array>0]);
				if(np.any(a_array<0)):
					x_max=np.min(temp[a_array<0]);
				
				while not x_max-x_min<1e-7*(np.abs(x_max)+np.abs(x_min)):
					x_new=0.5*(x_max+x_min);
					if (not np.isfinite(x_max)):
						x_new=2*np.abs(x_min+1);
					if (not np.isfinite(x_min)):
						x_new=-2*np.abs(x_max+1);
					Grad_x=(np.exp(x_new*a_array+b_array)-theta_array[:,i])*np.exp(x_new*a_array+b_array)*a_array;
					Grad_x=Grad_x.sum()
					if(Grad_x<0):
						x_min=x_new;
					else:
						x_max=x_new;
				
				
				
				Layer2[1,i]=0.5*(x_max+x_min);
				A_list[:,:,i]=Layer2[0,i]*I1+I2+Layer2[1,i]*I3
				
				conv[:,i]=np.dot(T,np.dot(A_list[:,:,i],Feature[:,i]));
	
				conv[:,i]=np.exp(conv[:,i]);
			
			
			
			#print "Layer2_1";
			
			LOSS=((Genome_BG[:,np.newaxis]*theta_array).sum(0)*P).sum()
			for i in range(4):
				temp=P[:,np.newaxis,np.newaxis,:]*theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:]
				temp=np.log(Genome_BG[np.newaxis,type[mid][i],np.newaxis]*temp.sum(3))
				LOSS-=(temp*X[:,type[mid][i],:]).sum()
			
			LOSS+=0.5*Reg*((theta_array-conv)**2).sum()
			
			#old1=new
			new=LOSS
			#if new>old1:
			#	print "wrong",new,old1;
			#	print Layer2;
			#	print old_layer;
			
			#print "finish";
		
		
		print r;

np.savetxt("C:\\MutationSignature\\comp\\base5_conv\\exp\\Feature_3.txt",Feature, fmt='%s',delimiter='\t')
np.savetxt("C:\\MutationSignature\\comp\\base5_conv\\exp\\Layer2_3.txt",Layer2, fmt='%s',delimiter='\t')
np.savetxt("C:\\MutationSignature\\comp\\base5_conv\\\exp\\Matrix_3.txt",Matrix, fmt='%s',delimiter='\t')
np.savetxt("C:\\MutationSignature\\comp\\base5_conv\\exp\\P_3.txt",P, fmt='%.18e',delimiter='\t')


conv=np.zeros((N,K));
	
for i in range(K):
	conv[:,i]=np.dot(T,np.dot(A_list[:,:,i],Feature[:,i]));

conv=np.exp(conv);
	
conv_all=np.zeros((N,3,K))

for i in range(4):
	conv_all[type[mid][i]]=theta_array[type[mid][i],np.newaxis,:]*Matrix[i,:,:];

conv_all=conv_all.reshape(N*3,K)
np.savetxt("C:\\MutationSignature\\comp\\base5_conv\\exp\\conv_3.txt",conv_all, fmt='%.18e',delimiter='\t')
