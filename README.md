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
