% algoritmo para treino da rede neural perceptron 1 neuronio e 1 camada
clc, clear;
load('data/classeA.mat');
load('data/classeB.mat');

w1 = rand();
w2 = rand();
w = [w1 w2];
teta = 1;
neta = 0.025;
total = [classeA classeB];
treino = [classeA(:,1:250) classeB(:,1:250)];
resposta_do_treino = [ones(1,250) ,-1*ones(1,250)];

erro =1;
epoca =0;
x = treino(1,:);
%reta incial antes do treino
reta_inicial = -(w1*x/w2)+(teta/w2);
%plot(reta_inicial,x,'r')
%hold on;

%init gif
figure(1)
filename = 'classificador.gif';
n=1;
%
while erro 
    n=n+1;
    erro=0;
    for i=1:500
        
        u1 = w1*treino(1,i);
        u2 = w2*treino(2,i);
        u = u1+u2-teta;
        if u>=0
            y=1;
        else
            y=-1;
        end
        if y~=resposta_do_treino(i)
           w1=w1 + neta*(resposta_do_treino(i)-y)*treino(1,i);
           w2=w2 + neta*(resposta_do_treino(i)-y)*treino(2,i);
           teta = teta +neta*(resposta_do_treino(i)-y)*-1;
           erro =1;
           
       % elseif y == resposta_do_treino(i)
        %    erro =0;    
        end
        
    end
    epoca = epoca+1;


    %plot(classeA(1,:),classeA(2,:),'ro')
    %hold on;
    %plot(classeB(1,:),classeB(2,:),'bo')
    %plot(treino(1,:),treino(2,:),'bx')
    %hold on;
      x = treino(1,:);
      reta = -(w1*x/w2)+(teta/w2);
      plot(reta,x,reta_inicial,x,'r', classeB(1,:),classeB(2,:),'bo',classeA(1,:),classeA(2,:),'ro', treino(1,:),treino(2,:),'bx')

      axis([-200 100 -200 100])
      
      %delay for gif 
      for j = 1:0.5:1000000
      end

      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end

%% plot result

      figure(2);
      reta = -(w1*x/w2)+(teta/w2);
      plot(reta,x,reta_inicial,x,'r', classeB(1,:),classeB(2,:),'bo',classeA(1,:),classeA(2,:),'ro', treino(1,:),treino(2,:),'bx')

      axis([-200 100 -200 100])

%% example for create a gif 
x = 0:0.01:1;
figure(3)
filename = 'create-gif.gif';
for n = 1:0.5:5
    
    for j = 1:0.5:1000000
    end
      y = x.^n;
      plot(x,y)
      drawnow
      frame = getframe(3);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end