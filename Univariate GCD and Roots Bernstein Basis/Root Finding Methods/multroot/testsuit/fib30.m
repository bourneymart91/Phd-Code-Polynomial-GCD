function [p,z] = fib30
%
%  test polynomial suggested by Goedecker
%
    n = 30;
    p = [-1,ones(1,n)];
    z = [-.9644254204388659, ...
        -.9439446810207413-.1984294062681023*i, ...
        -.9439446810207413+.1984294062681023*i, ...
        -.8833406981383851-.3885717512211873*i, ...
        -.8833406981383851+.3885717512211873*i, ...
        -.7850939960122262-.5624744046156477*i, ...
        -.7850939960122262+.5624744046156477*i, ...
        -.6532262984053521-.7128393846266873*i, ...
        -.6532262984053521+.7128393846266873*i, ...
        -.4931368697016981-.8333156264808382*i, ...
        -.4931368697016981+.8333156264808382*i, ...
        -.3113834310599106-.9187504784857605*i, ...
        -.3113834310599106+.9187504784857605*i, ...
        -.1154175582352257-.9653888511274815*i, ...
        -.1154175582352257+.9653888511274815*i, ...
        .8671260421202705e-1-.9710099019085114*i, ...
        .8671260421202705e-1+.9710099019085114*i, ...
        .2866757308842761-.9349927196966469*i, ...
        .2866757308842761+.9349927196966469*i, ...
        .4761639390711853-.8583045796885356*i, ...
        .4761639390711853+.8583045796885356*i, ...
        .6471568601471165-.7434108410885516*i, ...
        .6471568601471165+.7434108410885516*i, ...
        .7920805342279742-.5941279591518842*i, ...
        .7920805342279742+.5941279591518842*i, ...
        .9037637390940364-.4155315939426088*i, ...
        .9037637390940364+.4155315939426088*i, ...
        .9752028356220179-.2142999451221434*i, ...
        .9752028356220179+.2142999451221434*i, ...
        1.999999999068677];
    z = [z',ones(n,1)];
    if norm(imag(z(:,1))) == 0 
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.15f \t \t \t %3g \n', z');
    else
        fprintf('                 roots ')
        fprintf('   \t\t\t\t\t\t     multiplicities\n');
        fprintf('\n');
        fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
            [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;    