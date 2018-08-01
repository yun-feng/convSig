--main.lua

require "torch";

opt={};

dofile "/home/ted/conv/deep/multi/read_data_multisample.lua";
dofile "/home/ted/conv/deep/multi/model_multisample.lua";
--require "nn";
--model=torch.load("/home/ted/conv/deep/torch/adam/model2_bias_few_aver");

dofile "/home/ted/conv/deep/multi/criterion_multisample.lua";
dofile "/home/ted/conv/deep/multi/train_multisample.lua";

cycle=300000;

--test={}
--test.seq,test.sample,test.labels=readdata(45);
--old_loss=criterion:forward(model:forward(test.data),test.labels);

for c = 1,cycle do
        print(string.format("Start cycle %d", c));
        print("Loading data");
    train.seq,train.sample,train.labels=readdata(c%44);
        counter=0;
        print("Start train");
    training_loss=torch.Tensor(model_train());
        training_loss=training_loss:sum()/training_loss:size()[1]
        print(string.format("Average loss %6.6f",training_loss));
        --test_loss=criterion:forward(model:forward(test.data),test.labels);
       -- print(string.format("Test loss %6.6f",test_loss));
        print("Save model");
        if training_loss<old_loss then
             torch.save("/home/ted/conv/deep/multi/adam/model",model);
                old_loss=training_loss;
        end
end
