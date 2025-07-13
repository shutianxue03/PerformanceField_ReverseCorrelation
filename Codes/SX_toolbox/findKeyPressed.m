clc

CorrecKeyPressed = 0; pause(.5)
commandwindow
while ~CorrecKeyPressed
    [~, ~, keyCode] = KbCheck(-1); % gather all keyboards
    if sum(keyCode)>0,CorrecKeyPressed = 1; end
end

find(keyCode==1)