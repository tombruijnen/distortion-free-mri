function  res = PH (A,Id)
% Operator to split image into subimages with phase increments and the
% adjoint operator. 
% Forward: combine multiple images with phase adjustments
% Backward: Split images with phase adjustments according to A

res.A=A;
res.Id=Id;
res.adjoint=1; 
res=class(res,'PH');

%END
end