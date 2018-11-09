require "torch";

opt={};
train={}

dofile "/data/ted/multi/eso/deep/read_data_open_rc.lua";
dofile "/data/ted/multi/eso/deep/model_open_rc2.lua";


dofile "/data/ted/multi/eso/deep/criterion_open_rc.lua";
dofile "/data/ted/multi/eso/deep/train_open_rc.lua";

cycle=300000;


for c = 1,cycle do
        print(string.format("Start cycle %d", c));
        print("Loading data");
    train.seq,train.sample,train.labels,train.weights,train.mid=readdata(c%100);
        counter=0;
        print("Start train");
    training_loss=torch.Tensor(model_train());
        training_loss=training_loss:sum()/training_loss:size()[1]
        print(string.format("Average loss %6.6f",training_loss));
        print("Save model");
           torch.save("/data/ted//multi/eso/model_open_bg_rc2_400",model);
       
end
