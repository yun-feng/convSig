--feature extration

require "torch";
require "nn";
require "nngraph";

model=torch.load("/data/ted/multi/adam/model")

train={}

train.seq,train.sample,train.labels=readdata(0)

res=model:forward({train.seq,train.sample})

--model:forward(train.data)
sep=0.6

mask=torch.floor(torch.exp(model.output:select(2,2))/sep)

mask=mask+torch.floor(train.labels/2)
mask=torch.floor(mask/2)

idx=mask:nonzero()
idx=idx:view(idx:nElement())




signature=torch.Tensor({1,13,15,16,18,})


for s=1,signature:size(1) do
grad=torch.Tensor(train.labels:nElement(),20)
grad:zero()
grad:indexFill(2,torch.LongTensor{signature[s]},1)
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
out_file="./grad_"..s..".txt";
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


sep=0.36
mask=torch.floor(torch.exp(model.output:select(2,2))/sep)
mask=torch.ceil(mask/3)
mask=1-mask
mask=mask+1-torch.floor(train.labels/5)
mask=torch.floor(mask/2)

idx=mask:nonzero()
idx=idx:view(idx:nElement())

signature=torch.Tensor({2,3,4,5,6,7,9,10,11,14,17,19,20})

for s=1,signature:size(1) do
grad=torch.Tensor(train.labels:nElement(),20)
grad:zero()
grad:indexFill(2,torch.LongTensor{signature[s]},1)
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
out_file="./grad_low_"..s..".txt";
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



--save weights
save_tensor= model:get(14).weight

out = assert(io.open("./sig13_avh.txt", "w")) -- open a file for serialization

splitter = "\t"
for i=1,save_tensor:size(1) do
    for j=1,save_tensor:size(2) do
        out:write(save_tensor[i][j])
        if j == save_tensor:size(2) then
            out:write("\n")
        else
            out:write(splitter)
        end
    end
end

out:close()

wktensor=model:get(13).output
sig1=torch.ge((wktensor:t())[13],0.6)
sig1_neg=torch.le((wktensor:t())[13],0.2)
idx=sig1:nonzero()
idx=idx:view(idx:nElement())
sig1_av=train.seq:index(1,idx):sum(1)
sig1_av=sig1_av/idx:nElement()
save_tensor=sig1_av[1]

idx=sig1_neg:nonzero()
idx=idx:view(idx:nElement())
sig1_neg_av=train.seq:index(1,idx):sum(1)
sig1_neg_av=sig1_neg_av/idx:nElement()
save_tensor=sig1_neg_av[1]