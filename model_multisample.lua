--model.lua
require "torch";

require "nn";

require "nngraph";

Numbase=4;
SeqLen=1001; --from data.lua
Noutput=2;
Nsample=367

opt.NumKernels={8,12,16,20};


--mutation signature

---
---
--try to get the output size of the model
model=nn.Sequential();


model:add(nn.TemporalConvolution(Numbase, opt.NumKernels[1], 8 , 1)); --filter size: 8
model:add(nn.ELU(0.2));
model:add(nn.TemporalMaxPooling(4, 4)); --maxpooling 4 base



model:add(nn.TemporalConvolution(opt.NumKernels[1], opt.NumKernels[2], 8 , 1)); --filter size: 8
model:add(nn.ELU(0.2));
model:add(nn.TemporalMaxPooling(4, 4)); --maxpooling 4 base



model:add(nn.TemporalConvolution(opt.NumKernels[2], opt.NumKernels[3], 8 , 1)); --filter size: 8
model:add(nn.ELU(0.2));
model:add(nn.TemporalMaxPooling(4, 4));


model:forward(torch.rand(SeqLen,Numbase));
Ninput=model.output:size();
model:clearState();


model:forward(torch.rand(SeqLen,Numbase));
Ninput=model.output:size();
model:clearState();
---
---
---

--First CNN layer, CONV->ReLU->MaxPooling

h1= -nn.TemporalConvolution(Numbase, opt.NumKernels[1], 8 , 1); --filter size: 8
cnn_model= h1-nn.ELU(0.2);
cnn_model= cnn_model-nn.TemporalMaxPooling(4, 4); --maxpooling 4 base


--Second CNN layer
cnn_model= cnn_model-nn.TemporalConvolution(opt.NumKernels[1], opt.NumKernels[2], 8 , 1); --filter size: 8
cnn_model= cnn_model-nn.ELU(0.2);
cnn_model= cnn_model-nn.TemporalMaxPooling(4, 4); --maxpooling 4 base


--Third CNN layer
cnn_model= cnn_model-nn.TemporalConvolution(opt.NumKernels[2], opt.NumKernels[3], 8 , 1); --filter size: 8
cnn_model= cnn_model-nn.ELU(0.2);
cnn_model= cnn_model-nn.TemporalMaxPooling(4, 4);




--CNN Layer for Summary
cnn_model= cnn_model-nn.TemporalConvolution(opt.NumKernels[3], opt.NumKernels[4], Ninput[1] , 1); --sum to a single number
cnn_model= cnn_model-nn.ELU(0.2);
cnn_model= cnn_model-nn.AddConstant(0.2+1e-8);
cnn_model= cnn_model-nn.Reshape(opt.NumKernels[4]);

--Mutational Rate
h2= -nn.Linear(Nsample,5);
rate_model= h2-nn.Linear(5,opt.NumKernels[4]);
rate_model=rate_model -nn.Exp();


--Fully connected layer
full_model=nn.CMulTable()({cnn_model,rate_model})

full_model=full_model-nn.Replicate(Noutput,2);
full_model=full_model-nn.CMul(Noutput,opt.NumKernels[4]);


--Non-negative constraint
full_model=full_model-nn.ELU(0.2);
full_model=full_model-nn.AddConstant(0.2);

full_model=full_model-nn.Sum(3)

--Final normalize layer
full_model=full_model-nn.Normalize(1)
full_model=full_model-nn.Log()

model=nn.gModule({h1,h2}, {full_model})
