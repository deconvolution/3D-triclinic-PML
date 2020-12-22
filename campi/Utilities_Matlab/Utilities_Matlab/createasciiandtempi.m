clear all
listasac = textread('listasac.txt', '%s');
lista = textread('lista.txt', '%s');
tempi= zeros(length(listasac),4);
for i = 1:length(listasac)
    flag = 0;
    x1=0;
    x2=0;
    [t,data,SAChdr] = fget_sac(listasac{i});
    a = listasac{i};
    b1=a(36:47);
    b2=cat(2,a(12:17),'.',a(18:20));
    c1=str2double(b1(1:2));
    c2=str2double(b2(1:2));
    if c2-c1>0
        x1 = 60;
        flag = 300;
    end
    c3=str2double(b1(4:5));
    c4=str2double(b2(3:4));
    if c4-c3>0
        x2 = 60;
    end
    d1=(c4 - c3)*x2;
    c5=str2double(b1(7:12));
    c6=str2double(b2(5:9));
    t0=d1+c6-c5+flag;
    tempi(i,1)=t0;
    if SAChdr.times.a<0;
        tempi(i,2)=SAChdr.times.t0;
    else
        tempi(i,2)=SAChdr.times.a;
    end
    tempi(i,3)=SAChdr.times.t0;
    tempi(i,4)=SAChdr.times.t2;
    c=[t data];
    save(lista{i},'c','-ASCII')
end
for i = 1:length(listasac)
    if tempi(i,2)<0
        tempi(i,2)=300;
    end
end
for i = 1:length(listasac)
    if tempi(i,2)<0
        tempi(i,2)=300;
    end
    if tempi(i,1)>tempi(i,2) || isnan(tempi(i,1))||tempi(i,1)<0
        tempi(i,1) = tempi(i,2)-1;
    end
end
save -ascii tempi.txt tempi