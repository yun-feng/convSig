require "torch";

require "nn";

opt.weights={1,1,1}; --the training data is unbalanced for there are fewer hypermutated regions 

criterion=nn.ClassNLLCriterion(torch.Tensor(opt.weights));
