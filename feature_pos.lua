--feature extration

require "torch";
require "nn";
require "nngraph";

dofile "/data/ted/multi/eso/deep/read_data_open_rc.lua";

POS=tonumber(arg[1])
PROB=tonumber(arg[2])
PROB_L=torch.ceil(PROB*10)

sequence_file="/data/ted/multi/eso/channel_open_rc_si"..POS.."_"..PROB_L..".txt";
sample_file="/data/ted/multi/eso/sample_open_rc_si"..POS.."_"..PROB_L..".txt";
label_file="/data/ted/multi/eso/label_open_rc_si"..POS.."_"..PROB_L..".txt";
sequence_weight_file="/data/ted/multi/eso/weight_open.txt";
mid_file="/data/ted/multi/eso/mid_open_rc_si"..POS.."_"..PROB_L..".txt";


model=torch.load("/data/ted/multi/eso/model_open_bg_rc_si"..POS.."_"..PROB_L)

train={}

train.seq,train.sample,train.labels,train.weights,train.mid=readdata(0);

res=model:forward({train.seq,train.sample,train.mid})

--model:forward(train.data)


x=model:get(12).output
y=x:reshape(x:size(1),x:size(3))
--y=torch.exp(x)

Num_fil=10
for s=1,Num_fil do
sep=y:select(2,s):mean()+y:select(2,s):std()
mask=torch.floor(y:select(2,s)/sep)
idx=mask:nonzero()
idx=idx:view(idx:nElement())
grad=torch.Tensor(train.labels:nElement(),Num_fil)
grad:zero()
--ac
grad:select(2,s):fill(1)
--grad:indexCopy(2,torch.LongTensor{s},y:select(2,s):resize(y:size(1),1));
model:get(12):backward(model:get(11).output,grad);
for i = 2,11 do
	layer=12-i+1;
	model:get(layer):backward(model:get(layer-1).output,model:get(layer+1).gradInput);
	--if string.find(tostring(model.modules[layer]), 'ELU') then
	--    model:get(layer).gradInput=l:forward(model:get(layer).gradInput)
	--end
end
model:get(1):backward(train.seq,model:get(2).gradInput);
grad_t=model:get(1).gradInput:index(1,idx):sum(1)
grad_t=grad_t/idx:nElement()
seq_t=train.seq:index(1,idx):sum(1)/idx:nElement();
out_file="/data/ted/multi/eso/deep/bg_grad_rc_sm_si_"..s..POS.."_"..PROB_L..".txt";
save_tensor= grad_t
out = assert(io.open(out_file, "w")) -- open a file for serialization
splitter = "\t"
for i=1,save_tensor:size(2) do
for j=1,save_tensor:size(3) do
out:write(save_tensor[1][i][j])
if j == save_tensor:size(3) then
out:write("\n")
else
out:write(splitter)
end
end
end
out:close()
end

