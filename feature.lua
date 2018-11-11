--feature extration

require "torch";
require "nn";
require "nngraph";

dofile "/data/ted/multi/eso/deep/read_data_open_rc.lua";

model=torch.load("/data/ted/multi/eso/model_open_bg_rc2_400")

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
out_file="/data/ted/multi/eso/deep/bg_grad_rc2_sm_"..s..".txt";
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

save_tensor=model:get(19).weight
out = assert(io.open("/data/ted/multi/eso/deep/mean_mut_matrix.txt", "w")) -- open a file for serialization

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


--TCGA
--signature=torch.Tensor({1,13,15,16,18,})
--icgc
signature=torch.Tensor({4,7,8,9,10,11,12,14,15,16,18})

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
out_file="/data/ted/multi/coca/grad_"..s..".txt";
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

--signature=torch.Tensor({2,3,4,5,6,7,9,10,11,14,17,19,20})


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