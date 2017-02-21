function y=figsav(varargin)
%y=figsav(varargin). Сохранение графиков, varargin - список номеров графиков через запятую.
   for s=1:length(varargin)
      figure(varargin{s})
      if ~isempty(get(varargin{s},'Children'))
         rab2=0;
         while rab2==0
            Nff=input(['График ',num2str(varargin{s}),': введите имя файла сохранения графика (при отказе-<Enter>): '],'s');
            disp(' ')
            if (length(Nff)>1)
               if exist(Nff)==2 % m-file
                  disp('Файл с выбранным именем существует и будет уничтожен при')
                  disp('	сохранении  графиков с этим именем.')
                  rab3=input('Меняем имя файла сохранения? (у/<Enter>): ','s');
                  disp('')
                  if ~strcmp(rab3,'y') rab2=1; end
               else 
                  rab2=1;
               end
            else
               rab2=1;
            end
         end       %while
         v=axis;
         text(v(1),0.04*v(3)+0.96*v(4),[' ',Nff])
         eval(['print -f ',Nff,' -dmfile']);
      else
         close(varargin{s})
      end
   end
