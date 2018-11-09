--train.lua
require "torch"

require "math"

require "optim"


opt.KernelMax=0.9
opt.batch_size=60;
opt.epochs=1;

opt.State = {
   learningRate=0.001,
   learningRateDecay=1e-6,
   weightDecay=1e-7,
   beta1=0.9,
   beta2=0.99,
   epsilon=1e-8
}
opt.Method = optim.adam;

Num_sample=315;

model:training();

par,parGrad=model:getParameters();


feval=function(x)
    if x~=par then
        par:copy(x)
    end
    
    local start_index = counter * opt.batch_size + 1
    local end_index = math.min(seq_size[1]/4, (counter + 1) * opt.batch_size + 1)
    if end_index == seq_size[1]/4 then
        counter = 0
    else
        counter = counter + 1
    end
    local batch_inputs_seq = train.seq[{{start_index, end_index}, {}}]
	local batch_inputs_sample = train.sample[{{start_index, end_index}, {}}]
    local batch_targets = train.labels[{{start_index, end_index}}]
	local batch_weights = train.weights[{{start_index, end_index}}]
	local batch_mid = train.mid[{{start_index, end_index}}]    
    model:zeroGradParameters();
    
	--normalization for kernal 
    for i = 1,#model.modules do
        if string.find(tostring(model.modules[i]), 'TemporalConvolution') then
                model.modules[i].weight:renorm(2,1,opt.KernelMax)
        end
    end
	
	local p=torch.zeros(batch_inputs_seq:size(1),Num_sample);

		for i=1,batch_inputs_seq:size(1) do
                	p[i][(torch.floor(torch.rand(1)*Num_sample)+1)[1]]=1;
		end
                model:forward({batch_inputs_seq,p,batch_mid});
   local     f=0.045*torch.cmul(batch_weights,model.output[2]):sum()/model.output[2]:size(1);
                model:backward({batch_inputs_seq,p,bath_mid},{torch.zeros(batch_inputs_seq:size(1),Noutput),0.045*batch_weights/model.output[2]:size(1)});
	
    
     f=f+criterion:forward(model:forward({batch_inputs_seq,batch_inputs_sample,batch_mid})[1],batch_targets);
    
	f=f-1.0*(torch.log(model.output[2])):sum()/model.output[2]:size(1)	
    model:backward({batch_inputs_seq,batch_inputs_sample,bath_mid},{criterion:backward(model.output[1],batch_targets),-1.0/model.output[1]:size(1)*torch.cdiv(torch.ones(model.output[1]:size(1)),model.output[2])});
    
	
	
	
	
    return f,parGrad;
end


function model_train()
    
    local losses={}
    local iter=opt.epochs*math.ceil(seq_size[1]/4/opt.batch_size);
    
    for i = 1,iter do    
        local temp,BactchLoss=opt.Method(feval,par,opt.State);
        
        if i % 10 == 0 then -- don't print *every* iteration, this is enough to get the gist
            print(string.format("minibatches processed: %6s, loss = %6.6f", i, BactchLoss[1]))
        end
    
        losses[#losses + 1] = BactchLoss[1] -- append the new loss
    end
    
    return losses;
end
    
    
