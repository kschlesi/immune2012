       
    function [] = logsurf(data,Xaxis,Yaxis,format)
    
    % makes a surf plot with log-scale color map and colorbar
    % inputs: 'data' : m x n matrix of data to plot
    %         'Xaxis' : m x 1 vector giving x-axis values
    %         'Yaxis' : n x 1 vector givin y-axis values
    %         'format' (optional) : string to specify colorbar axis value
    %                               formatting; default = '%3.0G'
    
    % taking log of data (any neg. entries give -Inf)
    Dlog = log(data.*(data>0));
    
    % determining bounds of color axis labels
    Lex_max = floor(log10(max(max(data))));
    Lex_min = floor(log10(min(min(data(find(data>0)))))+1);
    L = transpose(10.^(Lex_min:Lex_max));
    
    % making figure
    figure
    surf(Xaxis,Yaxis,Dlog,'EdgeColor','none')
    axis([Xaxis(1) Xaxis(end) Yaxis(1) Yaxis(end)])
    logC = colorbar('Location','EastOutside');  
    formatSpec = '%3.0G';                       
    if nargin>3                                  
        formatSpec = format;
    end
    set(logC,'Ytick',log(L),'YTicklabel',num2str(L,formatSpec));
    
    end