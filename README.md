# methylation_nlp
Application of natural language processing techniques to predicting the methylome

20-12-21 
Run script window_func.py for cutting windows on methylation file from Tai.

15-02-22
Full implementation of lstm with connected neural net on bracewell. Major issue: data is mostly zeros, so model does very well by calling everything zero.
Solution: trying a categoriation 'zero or not'.
Speed also now an issue.

06-03-22
have now tried binary classification versus regression (aka single category). Have also tried 'balancing' the data - find all 1000bp windows with a methylated site. Match with random sample of 1000bp windows that do not contain a methylated site. Also tried only taking all 1000bp windows that do contain a methylated site. Still rapidly converges to zero. Data is also still zero inflated; 'balancing' in this context in non-trivial. 
Presented work to team. Comment was made that language has a natural grammar, whereas with DNA it's not obvious where the 'words' are. E.g. does it affect results if the boundaries between 10-mers are shifted by one base pair? I think this needs some serious thought, maybe language based models aren't as good a fit as we initially thought...

04-05-22
Will now try balancing by weighting a random selection of unmethylated windows to 1. All data is used to run the lstm but only methylated windows and a matching number of unmethylated windows used for the loss function.

10-05-22
Weighting via loss function now implemented. All evaluation metrics are terrible, but this is excellent because the model does not now converge to zero
instantly, and we can work on actually improving the model.
Probably need to give it more data at this point so will go back to working on scripts to format data. 
