disp 'hello'
filename = 'C:\Users\Home\Documents\NetBeansProjects\JavaFXApplication2\recorded.wav';
[y, fs] = wavread(filename);
y = y(:,1);
dt = 1/fs;
t = 0:dt:(length(y)*dt)-dt;
fid = fopen('C:\Users\Home\Documents\NetBeansProjects\JavaFXApplication2\amprecorded.txt','wt');  
fprintf(fid,'%d\n',y); 
fclose(fid);