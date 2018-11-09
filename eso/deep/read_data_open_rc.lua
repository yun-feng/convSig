require "torch";

sequence_file="/data/ted/multi/eso/channel_open_rc_si5.txt";
sample_file="/data/ted/multi/eso/sample_open_rc_si5.txt";
label_file="/data/ted/multi/eso/label_open_rc_si5.txt";
sequence_weight_file="/data/ted/multi/eso/weight_open400.txt";
mid_file="/data/ted/multi/eso/mid_open_rc_si5.txt";

Num_sam=315;

function readdata(x)
    
    local database = { }
    seqfile = io.open(sequence_file,"r")
    if seqfile then
        local counter=0;
        for seqline in seqfile:lines() do
                if((math.floor(counter/4))%100==x) then
                table.insert(database,seqline:split("\t"));
            end
        counter=counter+1;
        counter=counter%400;
            end
    end
      
	local mid={}
	midfile = io.open(mid_file,"r")
    if midfile then
        local counter=0;
        for midline in midfile:lines() do
                if((math.floor(counter))%100==x) then
                table.insert(mid,midline:split("\t"));
            end
        counter=counter+1;
        counter=counter%100;
            end
    end
	
	
	
        local labels={}
    labelfile= io.open(label_file,"r")
    if labelfile then
        local counter=0;
        for labelline in labelfile:lines() do
            if((math.floor(counter))%100==x) then
                table.insert(labels,labelline);
            end
        counter=counter+1;
        counter=counter%100;
        end
    end
        
    local samples = { }
        local sam_seq={};
    samfile = io.open(sample_file,"r");
    if samfile then
        local counter=0;
        for samline in samfile:lines() do
            if((math.floor(counter))%100==x) then
                sam_seq={}
		for i = 1,Num_sam do
                    sam_seq[i]=0
                end
                sam_seq[tonumber(samline)]=1;
                table.insert(samples,sam_seq);
            end
        counter=counter+1;
        counter=counter%100;
            end
    end
	
	local weights={};
	wtfile= io.open(sequence_weight_file,"r")
    if wtfile then
        local counter=0;
        for wtline in wtfile:lines() do
            if((math.floor(counter))%100==x) then
                table.insert(weights,wtline);
            end
        counter=counter+1;
        counter=counter%100;
        end
    end
	
		
	sequence=torch.Tensor(database);
    database=nil;
    seq_size=sequence:size();
    sequence:resize(seq_size[1]/4,4,seq_size[2]);
    sequence=sequence:transpose(2,3);	
	
    local shuffled_indices = torch.randperm(sequence:size(1)):long()
    sequence = sequence:index(1, shuffled_indices):squeeze()
    --length choosing
    --sequence=sequence:narrow(2,301,401);
    labels = torch.Tensor(labels):index(1, shuffled_indices):squeeze()
    mid = torch.Tensor(mid):index(1, shuffled_indices):squeeze()
    samples=torch.Tensor(samples):index(1, shuffled_indices):squeeze()
	weights = torch.Tensor(weights):index(1, shuffled_indices):squeeze()
        return sequence,samples,labels,weights,mid;
end
