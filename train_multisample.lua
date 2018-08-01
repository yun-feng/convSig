--train.lua
require "torch"

require "math"

require "optim"


opt.KernelMax=0.9
opt.batch_size=60;
opt.epochs=30;

opt.State = {
   learningRate=0.001,
   learningRateDecay=1e-5,
   weightDecay=1e-6,
   beta1=0.9,
   beta2=0.99,
   epsilon=1e-8
}
opt.Method = optim.adam;





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
    
    model:zeroGradParameters();
    
	--normalization for kernal 
    for i = 1,#model.modules do
        if string.find(tostring(model.modules[i]), 'TemporalConvolution') then
                model.modules[i].weight:renorm(2,1,opt.KernelMax)
        end
    end
    
    local f=criterion:forward(model:forward({batch_inputs_seq,batch_inputs_sample}),batch_targets);
    
    model:backward({batch_inputs_seq,batch_inputs_sample},criterion:backward(model.output,batch_targets));
    
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
    
    
