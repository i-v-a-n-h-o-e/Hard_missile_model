function y=figsav(varargin)
%y=figsav(varargin). ���������� ��������, varargin - ������ ������� �������� ����� �������.
   for s=1:length(varargin)
      figure(varargin{s})
      if ~isempty(get(varargin{s},'Children'))
         rab2=0;
         while rab2==0
            Nff=input(['������ ',num2str(varargin{s}),': ������� ��� ����� ���������� ������� (��� ������-<Enter>): '],'s');
            disp(' ')
            if (length(Nff)>1)
               if exist(Nff)==2 % m-file
                  disp('���� � ��������� ������ ���������� � ����� ��������� ���')
                  disp('	����������  �������� � ���� ������.')
                  rab3=input('������ ��� ����� ����������? (�/<Enter>): ','s');
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
