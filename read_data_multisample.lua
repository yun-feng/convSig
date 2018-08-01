--readdata.lua

require "torch";

sequence_file="/home/ted/conv/deep/multi/channel.txt";
sample_file="/home/ted/conv/deep/multi/TCGA_COAD_sample.txt";
label_file="/home/ted/conv/deep/multi/label.txt";



function readdata(x)
    --read in sequences
    local database = { }
    seqfile = io.open(sequence_file,"r")
    if seqfile then
        local counter=0;
        for seqline in seqfile:lines() do
	        if((math.floor(counter/4))%45==x) then
                table.insert(database,seqline:split("\t"));
            end
        counter=counter+1;
	    end
    end
	seqfile:close();
	--read in labels
	local labels={}
    labelfile= io.open(label_file,"r")
    if labelfile then
        local counter=0;
        for labelline in labelfile:lines() do
            if((math.floor(counter))%45==x) then
                table.insert(labels,labelline);
            end
        counter=counter+1;
        end
    end
	sequence=torch.Tensor(database);
	database=nil;
	seq_size=sequence:size();
    sequence:resize(seq_size[1]/4,4,seq_size[2]);
	sequence=sequence:transpose(2,3);
    --read in samples
    local samples = { }
    samfile = io.open(sample_file,"r")
    if samfile then
        local counter=0;
        for samline in samfile:lines() do
	        if((math.floor(counter))%45==x) then
                table.insert(samples,samline:split("\t"));
            end
        counter=counter+1;
	    end
    end
	local shuffled_indices = torch.randperm(sequence:size(1)):long()
    -- creates a shuffled *copy*, with a new storage
    sequence = sequence:index(1, shuffled_indices):squeeze()
    labels = torch.Tensor(labels):index(1, shuffled_indices):squeeze()
	labels=torch.floor(labels/4);
	samples=torch.Tensor(samples):index(1, shuffled_indices):squeeze()
	return sequence,samples,labels+1;
end


train={}





