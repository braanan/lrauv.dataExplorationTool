function [ result ] = isdeg( parameter )


if any( abs((pi/180) .* ( parameter )) > 1000 )
    result=true;
else
    result=false;
end

end