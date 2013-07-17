function a = flip(in)

  sizes = size(in);
  if (sizes(1) == 1)
    a = flipdim(in,2);
  elseif (sizes(2) == 1)
  a = flipdim(in,1);
  else
    error 'blah'
end
