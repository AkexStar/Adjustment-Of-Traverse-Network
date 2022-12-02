classdef daoxianwang_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ErrorAjustment      matlab.ui.Figure
        HomeMenu            matlab.ui.container.Menu
        MenuSave            matlab.ui.container.Menu
        Menu_2              matlab.ui.container.Menu
        Menu_8              matlab.ui.container.Menu
        Menu_9              matlab.ui.container.Menu
        Menu_10             matlab.ui.container.Menu
        BMenu               matlab.ui.container.Menu
        PMenu               matlab.ui.container.Menu
        lMenu               matlab.ui.container.Menu
        VMenu               matlab.ui.container.Menu
        xMenu               matlab.ui.container.Menu
        HelpMenu            matlab.ui.container.Menu
        Menu_4              matlab.ui.container.Menu
        Menu_5              matlab.ui.container.Menu
        AboutMenu           matlab.ui.container.Menu
        Menu_6              matlab.ui.container.Menu
        Menu_7              matlab.ui.container.Menu
        computeButton       matlab.ui.control.Button
        Label               matlab.ui.control.Label
        TextArea            matlab.ui.control.TextArea
        EditField_3Label    matlab.ui.control.Label
        A1EditField         matlab.ui.control.NumericEditField
        EditField_4Label    matlab.ui.control.Label
        A2EditField         matlab.ui.control.NumericEditField
        XmEditFieldLabel    matlab.ui.control.Label
        A1xEditField        matlab.ui.control.NumericEditField
        YmEditFieldLabel    matlab.ui.control.Label
        A1yEditField        matlab.ui.control.NumericEditField
        YmEditField_2Label  matlab.ui.control.Label
        A2yEditField        matlab.ui.control.NumericEditField
        XmEditField_2Label  matlab.ui.control.Label
        A2xEditField        matlab.ui.control.NumericEditField
        UIAxes              matlab.ui.control.UIAxes
        bianButton          matlab.ui.control.Button
        jiaoButton          matlab.ui.control.Button
        bianLamp            matlab.ui.control.Lamp
        jiaoLamp            matlab.ui.control.Lamp
        biantipEditField    matlab.ui.control.EditField
        jiaotipEditField    matlab.ui.control.EditField
        Label_2             matlab.ui.control.Label
        jiaosrcTextArea     matlab.ui.control.TextArea
        Label_3             matlab.ui.control.Label
        biansrcTextArea     matlab.ui.control.TextArea
        TabGroup            matlab.ui.container.TabGroup
        BTab                matlab.ui.container.Tab
        UITableB            matlab.ui.control.Table
        PTab                matlab.ui.container.Tab
        UITableP            matlab.ui.control.Table
        lTab                matlab.ui.container.Tab
        UITablel            matlab.ui.control.Table
        VTab                matlab.ui.container.Tab
        UITable             matlab.ui.control.Table
        xTab                matlab.ui.container.Tab
        UITablex            matlab.ui.control.Table
        Label_4             matlab.ui.control.Label
        Spinner             matlab.ui.control.Spinner
        Label_5             matlab.ui.control.Label
        EditFieldCount      matlab.ui.control.NumericEditField
    end

    
    properties (Access = public)
        bian % ÿÿÿ
        jiao % ÿÿÿ
        N_angl
        N_line
        N
        A1
        A2
        A1X
        A2X
        A1Y
        A2Y
        data
        count
    end
    
    methods (Access = public)
        
        function [fwj] = fwjjs(app,a1,a2,jiao,fwj,N_angl)
            for i=1:1:N_angl
                
                if (jiao(i,1)+1)==a1 && (jiao(i,2)+1)==a2 && fwj(a2,(jiao(i,3)+1))==0
                   	
                    fwj(a2,jiao(i,3)+1)=fwj(a1,a2)+jiao(i,4)+pi;
                    
                    fwj=fwjjs(app,a2,jiao(i,3)+1,jiao,fwj,N_angl);
                  		%ÿÿÿÿ
                  		
                elseif (jiao(i,3)+1)==a1 && (jiao(i,2)+1)==a2 && fwj(a2,(jiao(i,1)+1))==0
                   	%ÿÿÿÿÿÿa1 && ÿÿÿÿÿÿa2 && ÿÿÿÿa2,ÿÿÿÿÿÿÿÿÿ
                    fwj(a2,jiao(i,1)+1)=fwj(a1,a2)-jiao(i,4)+pi;
                    fwj=fwjjs(app,a2,jiao(i,1)+1,jiao,fwj,N_angl);
                  		
                elseif (jiao(i,1)+1)==a2 && (jiao(i,2)+1)==a1 && fwj(a1,(jiao(i,3)+1))==0
                   	%ÿÿÿÿÿÿa2 && ÿÿÿÿÿÿa1 && ÿÿÿÿa1,ÿÿÿÿÿÿÿÿÿ
                    fwj(a1,(jiao(i,3)+1))=fwj(a1,a2)+jiao(i,4);
                    fwj=fwjjs(app,a1,jiao(i,3)+1,jiao,fwj,N_angl);
                  		
                elseif (jiao(i,2)+1)==a1 && (jiao(i,3)+1)==a2 && fwj(a1,(jiao(i,1)+1))==0
                   	%ÿÿÿÿÿÿa2 && ÿÿÿÿÿÿa1 && ÿÿÿÿa1,ÿÿÿÿÿÿÿÿÿ
                    fwj(a1,(jiao(i,1)+1))=fwj(a1,a2)-jiao(i,4);
                    fwj=fwjjs(app,a1,(jiao(i,1)+1),jiao,fwj,N_angl);
                end
            end
        end
        function [X,Y] = zuobiao(app,A,AX,AY,bian,fwj,cc,X,Y,N_line)
           	X(A,1)=AX;
           	Y(A,1)=AY;
           	%ÿÿÿÿÿÿÿ
           	for i=1:1:N_line
              		if (bian(i,2)+1==A) && (cc(bian(i,3)+1,1)==0)
                 			cc(bian(i,3)+1,1)=1;
                 			X(bian(i,3)+1,1)=AX+bian(i,4)*cos(fwj(bian(i,2)+1,bian(i,3)+1));
                 			Y(bian(i,3)+1,1)=AY+bian(i,4)*sin(fwj(bian(i,2)+1,bian(i,3)+1));
                 			[X,Y]=zuobiao(app,bian(i,3)+1,X(bian(i,3)+1,1),Y(bian(i,3)+1,1),bian,fwj,cc,X,Y,N_line);
                 			
              		elseif(bian(i,3)+1==A) && (cc(bian(i,2)+1,1)==0)
                 			cc(bian(i,2)+1,1)=1;
                 			X(bian(i,2)+1,1)=AX+bian(i,4)*cos(fwj(bian(i,3)+1,bian(i,2)+1));
                 			Y(bian(i,2)+1,1)=AY+bian(i,4)*sin(fwj(bian(i,3)+1,bian(i,2)+1));
                 			[X,Y]=zuobiao(app,bian(i,2)+1,X(bian(i,2)+1,1),Y(bian(i,2)+1,1),bian,fwj,cc,X,Y,N_line);
              		end
           	end
        end
        function [B,bian,jiao,X01,v1,xxm,v,P,sgm0,sgmdianmm]= dxw(app,A1,A1X,A1Y,A2,A2X,A2Y,bian,jiao,N,N_angl,N_line)
           	%ÿÿÿÿÿ*****************************************
           	%A1ÿ	X	Y
           	%A2ÿ	X	Y
           	%ÿÿÿÿ
           	%bianÿ	ÿÿ	ÿÿÿ	ÿÿÿ	ÿÿÿÿ	ÿÿÿÿ
           	%jiaoÿ ÿÿÿ	ÿÿÿ	ÿÿÿ	ÿÿÿÿ	ÿÿÿÿ
           	%Nÿÿÿ
           	%N_anglÿÿÿÿÿ
           	%N_lineÿÿÿ
            
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿÿÿÿÿ
            %------------------------------------------------------------------------%
           	fwj=zeros(N,N);%ÿÿNÿÿ
           	%ÿÿÿ
           	%A1ÿ1ÿÿÿÿÿÿÿÿÿÿ0ÿÿÿmatlabÿÿÿÿÿ1ÿÿ
           	fwj(A1+1,A2+1)=atan2((A2Y-A1Y),(A2X-A1X));%ÿÿÿÿÿÿÿÿÿÿÿ
           	
           	fwj=fwjjs(app,A1+1,A2+1,jiao,fwj,N_angl);
           	for i=1:1:N
              		for j=1:1:N
                 			while fwj(i,j)>2*pi
                    				fwj(i,j)=fwj(i,j)-2*pi;
                 			end
                 			while    fwj(i,j)<0
                    				fwj(i,j)=fwj(i,j)+2*pi;
                 			end
              		end
           	end
            
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿÿÿ ÿÿÿÿÿ
            %------------------------------------------------------------------------%
           	%ÿÿ2ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ
           	X0=zeros(N,1);
           	Y0=zeros(N,1);
           	
           	%
           	X1=zeros(N,1);
           	Y1=zeros(N,1);
           	X2=zeros(N,1);
           	Y2=zeros(N,1);
           	%ÿÿÿ
           	global cc;
            cc=zeros(N,1);
           	[X1,Y1]=zuobiao(app,A2+1,A2X,A2Y,bian,fwj,cc,X1,Y1,N_line);
           	%ÿÿÿÿÿ*****************************************
           	%ÿÿÿÿÿ	ÿÿÿX	ÿÿÿY ÿÿÿÿÿ	ÿÿÿÿÿ	ccÿÿÿÿ	ÿÿXÿÿÿ	ÿÿYÿÿÿ
           	cc=zeros(N,1);
           	[X2,Y2]=zuobiao(app,A1+1,A1X,A1Y,bian,fwj,cc,X2,Y2,N_line);
           	
           	for i=1:1:N
               	%ÿÿNÿÿ	ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ ÿÿÿÿ
              		X0(i,1)=(X1(i,1)+X2(i,1))/2;
              		Y0(i,1)=(Y1(i,1)+Y2(i,1))/2;
           	end
           	X0(A1+1)=A1X;
           	Y0(A1+1)=A1Y;
           	X0(A2+1)=A2X;
           	Y0(A2+1)=A2Y;
           	%ÿÿÿÿÿÿÿ
            
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿÿÿÿÿÿÿ Nÿÿ ÿÿÿÿX,Y ÿÿN*2
            %------------------------------------------------------------------------%
           	x0=zeros(N*2,1);
           	
           	for i=1:1:N
              		x0(2*i-1,1)=X0(i,1);%ÿÿÿÿX
              		x0(2*i,1)=Y0(i,1);	%ÿÿÿÿY
           	end
            
           	x0([2*A1+1,2*A1+2,2*A2+1,2*A2+2],:)=[];%ÿÿÿÿÿÿ
           	%ÿÿÿÿÿÿX,Yÿÿÿÿÿ
            
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿÿÿÿÿÿÿÿÿ jiaoxishu
            %------------------------------------------------------------------------%
           	%ÿÿN_anglÿÿÿÿ
           	jiaoxishu=zeros(N_angl,N*2);
           	%ÿÿjiaoxishu
           	
           	for i=1:1:N_angl
              		j=jiao(i,2)+1;%ÿÿÿ ÿÿj
              		k=jiao(i,3)+1;%ÿÿÿ ÿÿk
              		h=jiao(i,1)+1;%ÿÿÿ ÿÿh
              		
              		rou=180*3600/pi;
              		%ÿÿÿÿÿÿÿÿÿ
              		
              		up1=rou*sin(fwj(j,k))/1000;
              		up2=-rou*cos(fwj(j,k))/1000;
              		%ÿÿÿ ÿ/mm ÿÿÿÿ1000
              		
              		under1=sqrt((X0(k,1)-X0(j,1))^2+(Y0(k,1)-Y0(j,1))^2);
              		%S0 ÿÿÿÿ
              		
              		a1=up1/under1;
              		b1=up2/under1;
              		%ajk ÿ bjk
              		
              		up3=rou*sin(fwj(j,h))/1000;
              		up4=-rou*cos(fwj(j,h))/1000;
              		
              		under2=sqrt((X0(h,1)-X0(j,1))^2+(Y0(h,1)-Y0(j,1))^2);
              		
              		a2=up3/under2;
              		b2=up4/under2;
              		
              		%iÿÿÿÿÿ
              		
              		jiaoxishu(i,2*h-1)=a2;	%hÿxÿÿ
              		jiaoxishu(i,2*h)=b2;	%hÿYÿÿ
              		
              		jiaoxishu(i,2*j-1)=a1-a2;
              		jiaoxishu(i,2*j)=b1-b2;
              		
              		jiaoxishu(i,2*k-1)=-a1;
              		jiaoxishu(i,2*k)=-b1;
              		
           	end
           	
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿÿÿÿ bianxishu
            %------------------------------------------------------------------------%
           	%ÿÿN_lineÿÿ
           	bian0=zeros(N_line,1);%ÿÿÿÿÿÿÿbian0
           	for  i=1:1:N_line
              		j1=bian(i,2);
              		%ÿÿÿ
              		j2=bian(i,3);
              		%ÿÿÿ
              		bian0(i)=sqrt((X0(j1+1)-X0(j2+1))*(X0(j1+1)-X0(j2+1))+(Y0(j1+1)-Y0(j2+1))*(Y0(j1+1)-Y0(j2+1)));
           	end
           	
           	bianxishu=zeros(N_line,N*2);%ÿÿÿÿÿÿÿÿÿÿ
           	for i=1:1:N_line
              		j=bian(i,2);%ÿÿÿÿÿ
              		k=bian(i,3);%ÿÿÿÿÿ
              		up1=X0(k+1)-X0(j+1);	%ÿXjk
              		up2=Y0(k+1)-Y0(j+1);	%ÿYjk
              		
              		under=bian0(i);	%ÿÿÿÿÿ
              		m=up1/under;	%ÿX
              		n=up2/under;	%ÿY
              		bianxishu(i,2*(j+1)-1)=-m;
              		bianxishu(i,2*(j+1))=-n;
              		bianxishu(i,2*(k+1)-1)=m;
              		bianxishu(i,2*(k+1))=n;
           	end
           	
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿÿÿB
            %------------------------------------------------------------------------%
           	B=[jiaoxishu;bianxishu];
           	B(N_angl+N_line,:)=[];
           	%ÿÿÿN_angl+N_lingÿÿÿÿÿ
           	B(:,[2*A1+1,2*A1+2,2*A2+1,2*A2+2])=[];
           	%ÿÿÿÿÿÿÿÿÿÿÿÿÿ
           	
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿPÿÿ
            %------------------------------------------------------------------------%
           	P=zeros(N_angl+N_line,N_angl+N_line);
           	%N_angl+N_line ÿÿÿÿ
           	for i=1:1:N_angl
              		P(i,i)=1;
              		%ÿÿÿÿÿÿÿ1
           	end
           	j=1;
           	for i=N_angl+1:1:N_angl+N_line
              		P(i,i)=4*12*12/bian(j,5)/bian(j,5);
              		%ÿÿÿÿÿ
              		j=j+1;
           	end
            
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿlÿÿ
            %------------------------------------------------------------------------%
           	l1=zeros(N_line,1);
           	%ÿÿÿÿl1ÿÿÿ
           	for i=1:N_line
              		%bian0 ÿÿÿÿ
              		%bian(i,5) ÿÿÿÿ
              		l1(i,1)=(bian(i,5)-bian0(i))*1000;%ÿÿÿÿÿÿ
              		if l1(i,1)<0.0000001 && l1(i,1)>-0.0000001
                 			l1(i,1)=0;
                 			%ÿÿÿÿÿÿ ÿÿ0
              		end
           	end
           	
           	l1(N_line,:)=[];
           	%ÿÿÿÿÿÿÿÿ
           	
           	jiao01=zeros(N_angl,1);
           	jiao02=zeros(N_angl,1);
           	jiao0=zeros(N_angl,1);
           	
           	for i=1:1:N_angl %ÿÿÿÿÿÿÿÿÿÿÿÿÿ
              		jiao01(i,1)=atan2((Y0(jiao(i,3)+1,1)-Y0(jiao(i,2)+1,1)),(X0(jiao(i,3)+1,1)-X0(jiao(i,2)+1,1)));
              		%jkÿÿÿ
              		jiao02(i,1)=atan2((Y0(jiao(i,1)+1,1)-Y0(jiao(i,2)+1,1)),(X0(jiao(i,1)+1,1)-X0(jiao(i,2)+1,1)));
              		%jhÿÿÿ
              		%atan2 ÿÿ Y ÿ X ÿÿÿÿÿÿÿ
              		jiao0(i,1)=jiao01(i,1)-jiao02(i,1);
              		%jk-jh=hjk
           	end
           	
           	for i=1:1:N_angl %ÿÿÿÿÿÿ 0-2*pi
                while jiao0(i,1)>2*pi
                    jiao0(i,1)=jiao0(i,1)-2*pi;
                end
                while jiao0(i,1)<0
                    jiao0(i,1)=jiao0(i,1)+2*pi;
                end
           	end
           	
           	
           	l2=zeros(N_angl,1);
           	for i=1:1:N_angl     %ÿÿÿÿÿÿÿÿlÿÿÿÿÿÿÿÿ
              		l2(i,1)=(jiao(i,5)-jiao0(i))*3600*180/pi;
              		if l2(i,1)<0.00000001 && l2(i,1)>-0.00000001
                 			l2(i,1)=0;
              		end
           	end
           	
           	l=[l2;l1];%ÿÿÿÿl
           	
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿÿÿ
            %------------------------------------------------------------------------%
           	P(N_angl+N_line,:)=[];
           	P(:,N_angl+N_line)=[];
           	%ÿÿÿÿÿÿÿÿÿÿ
           	NBB=B'*P*B;
           	W=B'*P*l;
           	xx=NBB\W;%ÿÿmm
           	v=B*xx-l;
           	%ÿÿÿÿÿÿv ÿÿÿmm
           	
           	xxm=xx.*0.001;
           	%ÿÿÿÿÿÿÿÿÿÿÿ m
            
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿ
            %------------------------------------------------------------------------%
           	sgm0=sqrt((v'*P*v)/(N_angl+N_line-1-(N*2-4)));
           	%ÿÿÿÿÿÿ n-t=n-u
           	X01=x0+xxm;%ÿÿÿÿÿ
           	
           	v1=zeros(N_line+N_angl,1);%ÿÿÿÿ
           	for i=N_angl+1:1:N_angl+N_line-1
              		v1(i,1)=bian(i-N_angl,5)+v(i,1)*0.001;
           	end
           	v1(N_angl+N_line,1)=bian(N_line,5);
           	
           	for i=1:1:N_angl
              		v1(i,1)=jiao(i,5)+v(i,1)/3600*pi/180;
           	end
           	
           	Qxx=inv(NBB);%ÿÿÿÿ
           	Qx=zeros(N-2,1);%ÿxÿÿÿÿÿÿ
           	for i=1:1:N-2
              		Qx(i,1)=Qxx(2*i-1,2*i-1);
           	end
           	Qy=zeros(N-2,1);%ÿyÿÿÿÿÿÿ
           	for i=1:1:N-2
              		Qy(i,1)=Qxx(2*i,2*i);
           	end
           	
           	sgmdian=zeros(N-2,1);
           	%ÿÿÿÿÿÿÿÿ ÿÿÿmm
           	for i=1:1:N-2
              		sgmdian(i,1)=sgm0*sqrt(Qx(i,1)+Qy(i,1));
           	end
           	
           	sgmdianmm=zeros(N-2,1);
           	%ÿÿÿÿÿÿÿÿ ÿÿÿm
           	for i=1:1:N-2
              		sgmdianmm(i,1)=sgmdian(i,1);
           	end
           	
            %************************************************************************%
            %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
           	%ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ
            %------------------------------------------------------------------------%
           	for i=1:1:N_line-1
              		bian(i,4)=bian(i,5)+v(i+N_angl,1)*0.001;
           	end
           	for i=1:1:N_angl
              		jiao(i,4)=v(i,1)/3600*pi/180+jiao(i,5);
            end
            
            
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.count=0;
            app.bianLamp.Color="red";
                app.jiaoLamp.Color="red";
                app.biantipEditField.Value="ÿÿÿÿÿÿÿ";
                app.jiaotipEditField.Value="ÿÿÿÿÿÿÿ";
                app.A1EditField.Value=0;
                app.A2EditField.Value=0;
                app.A1xEditField.Value=0;
                app.A2xEditField.Value=0;
                app.A1yEditField.Value=0;
                app.A2yEditField.Value=0;
                app.N=0;
                app.bian=[];
                app.jiao=[];
                app.N_angl=0;
                app.N_line=0;
        end

        % Button pushed function: computeButton
        function computeButtonPushed(app, event)
            if app.bianLamp.Color(1,1)==1 || app.jiaoLamp.Color(1,1)==1
                app.TextArea.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";
                return;
            end
            %ÿÿÿÿÿÿÿ ÿÿÿÿÿÿÿÿÿ
            try
                app.A1=app.A1EditField.Value;%A1ÿÿ
                app.A2=app.A2EditField.Value;%A2ÿÿ
                %A1 A2ÿÿ
                app.A1X=app.A1xEditField.Value;
                app.A1Y=app.A1yEditField.Value;
                app.A2X=app.A2xEditField.Value;
                app.A2Y=app.A2yEditField.Value;
                
                app.N=max(app.jiao(:,1:3),[],'all','omitnan')+1;
                %Nÿÿÿ
                app.N_angl=length(app.jiao(:,1));
                %ÿÿÿÿÿ
                
                knownLine=sqrt((app.A1X-app.A2X)^2+(app.A1Y-app.A2Y)^2);
                app.N_line=length(app.bian(:,1));
                app.bian=[[1:1:app.N_line]',app.bian];
                %ÿÿÿÿÿ ÿÿÿÿÿ
                app.bian=[app.bian;[app.N_line+1, app.A1, app.A2, knownLine, knownLine]];
                %ÿÿÿÿÿ
                app.N_line=app.N_line+1;
                %ÿÿÿ
                
                i=1;
                app.count=1;
                while i==1%ÿÿÿÿÿÿ
                   	%dxw ÿÿÿÿÿÿÿ
                   	[~,app.bian,app.jiao,~,~,xxm1,~,~,~,~]=dxw(app,app.A1,app.A1X,app.A1Y,app.A2,app.A2X,app.A2Y,app.bian,app.jiao,app.N,app.N_angl,app.N_line);
                    
                   	%ÿÿÿÿÿ
                   	%A1ÿ	X	Y
                   	%A2ÿ	X	Y
                   	%bianÿÿ	ÿÿ	ÿÿÿ	ÿÿÿ	ÿÿÿÿÿÿ	ÿÿÿÿÿÿ
                   	%jiaoÿÿ ÿÿÿ	ÿÿÿ	ÿÿÿ	ÿÿÿÿÿÿ	ÿÿÿÿÿÿ
                    %N ÿÿÿ
                    %N_angl ÿÿÿÿÿ
                    %N_line ÿÿÿ
                    
                   	%ÿÿÿÿÿ
                   	%bian
                   	%jiao
                   	%Xÿÿÿÿÿ
                   	%L	ÿÿÿÿÿÿÿÿÿ
                   	%xxm ÿÿÿÿÿ ÿmÿÿÿ
                   	%V	ÿÿÿ
                   	%P	ÿÿ
                   	%sgm0 ÿÿÿÿÿÿ
                    %sgmdianmi	ÿÿÿÿÿÿÿÿ ÿÿÿm
                    
                   	[B,app.bian,app.jiao,X,L,xxm2,V,P,sgm0,sgmdianmm]=dxw(app,app.A1,app.A1X,app.A1Y,app.A2,app.A2X,app.A2Y,app.bian,app.jiao,app.N,app.N_angl,app.N_line);
                    app.data{app.count,1}=B;
                    app.data{app.count,2}=P;
                    l=B*xxm2-V;
                    app.data{app.count,3}=l;
                    app.data{app.count,4}=V;
                    app.data{app.count,5}=xxm2;
                    app.count=app.count+1;
                    if max(abs(xxm2-xxm1))>0.00000001
                      		i=1;
                    else
                      		i=0;
                    end
                end
                
                Nbb=B'*P*B;
                Qll=B*inv(Nbb)*B';
                Dll=sgm0^2*Qll;
                Dll2=diag(Dll);
                sgmll=sqrt(Dll2);
                [sgmllmax,index]=max(sgmll(app.N_angl+1:app.N_angl+app.N_line-1,1));
                %ÿÿÿÿÿÿÿ ÿÿÿmm
                sgmllmax=sgmllmax/1000/app.bian(index,4);
                %ÿÿÿÿÿÿÿÿÿ
                
                L(1:app.N_angl,1)=rad2deg(L(1:app.N_angl,1));
                %ÿÿÿÿÿÿÿÿÿÿÿÿ
                dushu1=L(1:app.N_angl,1);
                dushu2=string(fix(dushu1))+"°";
                dushu1=(dushu1-fix(dushu1))*60;
                dushu2=dushu2+string(fix(dushu1))+"ÿ";
                dushu1=(dushu1-fix(dushu1))*60;
                dushu1=round(dushu1,1);
                %ÿÿÿÿÿÿÿÿÿÿ1ÿ
                dushu(:,2)=dushu2+string(dushu1)+"ÿ";
                dushu(:,1)="ÿ"+string(app.jiao(:,1))+string(app.jiao(:,2))+string(app.jiao(:,3))+"=";
                %dushuÿÿÿÿÿÿÿ ÿÿÿÿÿÿ
                
                L(app.N_angl+1:app.N_angl+app.N_line,1)=round(L(app.N_angl+1:app.N_angl+app.N_line,1),3);
                %L(N_angl+1:N_angl+N_line,1)ÿÿÿÿÿÿÿÿÿÿÿÿÿ3ÿÿÿÿÿm
                
                V(1:app.N_angl)=round(V(1:app.N_angl),1);
                dushu(:,3)=string(V(1:app.N_angl,1));
                %dushu(:,3)ÿÿÿÿÿÿÿÿÿÿÿÿ1ÿÿÿÿÿÿ
                V(1+app.N_angl:app.N_angl+app.N_line-1,1)=round(V(1+app.N_angl:app.N_angl+app.N_line-1,1)/1000,3);
                %ÿÿÿÿÿÿÿÿÿÿÿÿ3ÿÿÿÿÿm
                bianchang(:,1)=string(app.bian(:,2));
                bianchang(:,2)=string(app.bian(:,3));
                bianchang(:,3)=string(L(app.N_angl+1:app.N_angl+app.N_line,1));
                bianchang(1:app.N_line-1,4)=string(V(app.N_angl+1:app.N_angl+app.N_line-1,1));
                %bianchangÿÿÿÿÿÿÿÿÿÿ ÿÿÿ
                
                X=round(X,3);
                %ÿÿÿÿÿÿÿÿÿÿÿ
                tempX=zeros(app.N-2,1);
                tempY=zeros(app.N-2,1);
                cla(app.UIAxes);
                %ÿÿÿÿÿ
                hold (app.UIAxes,"on");
                for i=1:app.N-2
                    tempX(i)=round(X(2*i-1,1),3);
                    tempY(i)=round(X(2*i,1),3);
                end
                coor=zeros(app.N,2);
                coor(app.A1+1,1)=app.A1X;
                coor(app.A1+1,2)=app.A1Y;
                coor(app.A2+1,1)=app.A2X;
                coor(app.A2+1,2)=app.A2Y;
                j=1;
                for i=0:app.N-3
                    if j-1==app.A1
                        j=j+1;
                        %ÿÿÿÿÿ
                    end
                    if j-1==app.A2
                        j=j+1;
                        %ÿÿÿÿÿ
                    end
                    coor(j,1)=X((i+1)*2-1,1);
                    coor(j,2)=X((i+1)*2,1);
                    j=j+1;
                end
                %coor ÿÿÿÿÿÿÿÿ ÿÿÿÿÿÿÿÿÿÿ
                plot(app.UIAxes,tempY,tempX,'bo');
                plot(app.UIAxes,app.A1Y,app.A1X,'r^');
                plot(app.UIAxes,app.A2Y,app.A2X,'r^');
                %ÿÿ
                for i=1:app.N_line
                    plot(app.UIAxes,[coor(app.bian(i,2)+1,2),coor(app.bian(i,3)+1,2)],[coor(app.bian(i,2)+1,1),coor(app.bian(i,3)+1,1)],'-c');
                end
                %ÿÿ
                for i=1:app.N
                    text(app.UIAxes,coor(i,2)+5,coor(i,1)+5,string(i-1));
                end
                %ÿÿÿÿ
                
                %ÿÿÿÿÿÿÿ
                app.TextArea.Value="*************ÿÿÿÿÿ*************";
                for i=1:app.N_angl
                    app.TextArea.Value=[app.TextArea.Value;dushu(i,1)+dushu(i,2)];
                end
                app.TextArea.Value=[app.TextArea.Value;"*************ÿÿÿÿÿ*************"];
                for i=1:app.N_angl
                    app.TextArea.Value=[app.TextArea.Value;dushu(i,3)];
                end
                app.TextArea.Value=[app.TextArea.Value;"*************ÿÿÿÿÿ*************"];
                for i=1:app.N_line-1
                    app.TextArea.Value=[app.TextArea.Value;bianchang(i,1)+"--"+bianchang(i,2)+"    "+bianchang(i,3)];
                end
                app.TextArea.Value=[app.TextArea.Value;"*************ÿÿÿÿÿ*************"];
                for i=1:app.N_line-1
                    app.TextArea.Value=[app.TextArea.Value;bianchang(i,1)+"--"+bianchang(i,2)+"    "+bianchang(i,4)];
                end
                app.TextArea.Value=[app.TextArea.Value;"*************ÿÿÿÿÿ*************"];
                for i=1:(app.N)
                    app.TextArea.Value=[app.TextArea.Value;string(i-1)+"  X:"+string(coor(i,1))+"  Y:"+string(coor(i,2))];
                end
                app.TextArea.Value=[app.TextArea.Value;"*************ÿÿÿÿÿ*************"];
                for i=1:(app.N-2)
                    app.TextArea.Value=[app.TextArea.Value;string(round(sgmdianmm(i,1),3))];
                end
                
                app.TextArea.Value=[app.TextArea.Value;"*************ÿÿÿÿÿÿÿÿÿÿÿ*************"];
                app.TextArea.Value=[app.TextArea.Value;string(round(sgmllmax,5))+" (1/"+string(round(1/sgmllmax,0))+")"];
                app.TextArea.Value=[app.TextArea.Value;"*************ÿÿÿÿÿÿÿÿÿÿÿ*************"];
                app.TextArea.Value=[app.TextArea.Value;string(round(max(sgmdianmm,[],'all','omitnan'),3))];
                app.TextArea.Value=[app.TextArea.Value;"*************ÿÿÿÿÿÿÿÿÿ*************"];
                app.TextArea.Value=[app.TextArea.Value;string(round(sgm0,3))];
                
                app.Spinner.Limits=[1,app.count-1];
                app.EditFieldCount.Value=app.count-1;
                value=1;
                app.UITableB.Data=app.data{value,1};
                app.UITableP.Data=app.data{value,2};
                app.UITablel.Data=app.data{value,3};
                app.UITable.Data=app.data{value,4};
                app.UITablex.Data=app.data{value,5};
                
                %ÿÿÿÿÿÿÿ ÿÿÿÿÿÿÿÿ
                app.bianLamp.Color="red";
                app.jiaoLamp.Color="red";
                app.biantipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";
                app.jiaotipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";
                app.A1EditField.Value=0;
                app.A2EditField.Value=0;
                app.A1xEditField.Value=0;
                app.A2xEditField.Value=0;
                app.A1yEditField.Value=0;
                app.A2yEditField.Value=0;
                app.N=0;
                app.N_angl=0;
                app.N_line=0;
                app.bian=[];
                app.jiao=[];
            catch ME
                app.bianLamp.Color="red";
                app.jiaoLamp.Color="red";
                app.biantipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";
                app.jiaotipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";
                app.A1EditField.Value=0;
                app.A2EditField.Value=0;
                app.A1xEditField.Value=0;
                app.A2xEditField.Value=0;
                app.A1yEditField.Value=0;
                app.A2yEditField.Value=0;
                app.N=0;
                app.bian=[];
                app.jiao=[];
                app.N_angl=0;
                app.N_line=0;
                throw(ME);
            end
            %app.UITableB.Data=dushu;
            %app.TextArea.Value = uitextarea('Value',temp);
            %             app.TextArea.Value=app.TextArea.Value+mat2str(dushu);
            %             app.TextArea.Value=app.TextArea.Value+"\n";
            %             app.TextArea.Value=app.TextArea.Value+mat2str(bianchang);
        end

        % Callback function
        function bianEditFieldValueChanging(app, event)
            
        end

        % Callback function
        function jiaoEditFieldValueChanging(app, event)
            
        end

        % Menu selected function: MenuSave
        function MenuSaveSelected(app, event)
            A = app.TextArea.Value;
            [file,path] = uiputfile('*.xlsx');
            if path==0
                return;
            end
            filename = fullfile(path,file);
            xlswrite(filename,A);
            
        end

        % Menu selected function: Menu_2
        function Menu_2Selected(app, event)
            h=get(app.UIAxes,'children');
            figure('visible','off');
            copyobj(h,gca);
            [file,path] = uiputfile('*.png');
            if path==0
                return
            end
            filename = fullfile(path,file);
            saveas(gcf,filename);
        end

        % Menu selected function: Menu_6
        function Menu_6Selected(app, event)
            app.TextArea.Value=["ÿÿÿÿÿÿÿÿ ©ÿÿÿ";"ÿÿÿÿÿÿÿÿÿÿÿÿÿ";"ÿÿÿÿÿ";"ÿÿÿÿÿÿÿÿÿ";"ÿÿÿÿÿÿÿÿÿÿÿÿ";"©Li Jintao, School of Geodesy and Geomatics, Wuhan University"];
        end

        % Menu selected function: Menu_7
        function Menu_7Selected(app, event)
            app.TextArea.Value=["Email:";"lijintao@whu.edu.cn";"Github:";"https://github.com/Alex-Starry-2019"];
        end

        % Menu selected function: Menu_5
        function Menu_5Selected(app, event)
            app.TextArea.Value=["ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿlijintao@whu.edu.cn";"ÿÿÿÿÿÿÿÿÿÿ"];
        end

        % Menu selected function: Menu_4
        function Menu_4Selected(app, event)
            app.TextArea.Value=[
                "ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿV2.0";...
                "ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";...
                "ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";...
                "";...
                "ÿÿÿÿÿÿÿ";...
                "ÿ1ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";...
                "ÿ2ÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";...
                "ÿ3ÿÿÿÿÿÿÿÿ ÿÿÿÿÿÿÿÿÿÿ2ÿÿÿÿÿÿ";...
                "ÿ4ÿÿÿÿÿÿÿÿÿÿ";...
                "ÿ5ÿÿÿÿÿÿÿÿÿ";...
                "";...
                "ÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";...
                "ÿÿÿÿÿÿÿÿÿÿÿÿ(m),ÿÿÿÿÿÿÿÿÿÿÿÿ";...
                "0 -- 1	97.711";...
                "1 -- 2	69.81";...
                "2 -- 3	78.536";...
                "3 -- 4	72.12";...
                "4 -- 5	67.714";...
                "5 -- 6	ÿÿÿ";...
                "6 -- 0	65.505";...
                "3 -- 7	77.436";...
                "7 -- 8	40.1";...
                "8 -- 6	68.199";...
                "";...
                "ÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";...
                "ÿÿÿÿÿÿÿÿÿÿÿ(°ÿÿ),ÿÿÿÿÿÿÿÿÿÿÿ:";...
                "ÿ0 1 2	236°00ÿ33.5ÿ";...
                "ÿ1 0 6	130°33ÿ18.9ÿ";...
                "ÿ1 2 3	225°51ÿ8.5ÿ";...
                "ÿ2 3 4	231°45ÿ21ÿ";...
                "ÿ2 3 7	281°25ÿ11.9ÿ";...
                "ÿ3 4 5	213°46ÿ6.1ÿ";...
                "ÿ4 3 7	49°40ÿ8.4ÿ";...
                "ÿ4 5 6	246°24ÿ48ÿ";...
                "ÿ5 6 0	236°45ÿ44.4ÿ";...
                "ÿ5 6 8	323°59ÿ33.9ÿ";...
                "ÿ0 6 8	87°13ÿ30.7ÿ";...
                "ÿ3 7 8	174°53ÿ58.8ÿ";...
                "ÿ7 8 6	199°35ÿ58.9ÿ";...
                "";...
                "ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";...
                "5 x=164.668(m), y=112.313(m)";...
                "6 x=274.722(m), y=136.706(m)";...
                "";...
                "ÿÿÿÿÿÿÿÿÿÿ1/2000";...
                "ÿÿÿÿÿÿÿ12ÿ";...
                ];
        end

        % Callback function
        function SwitchValueChanged(app, event)
            
        end

        % Button pushed function: bianButton
        function bianButtonPushed(app, event)
            
            if app.biansrcTextArea.Value==""
                app.bianLamp.Color="red";
                app.biantipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";
                return;
            end
            
            value=strings;
            for i=1:length(app.biansrcTextArea.Value)
                value=value+" "+app.biansrcTextArea.Value{i,1};
            end
            bian_temp1 = textscan(value,'%f %s %f %s');
            [~,n]=size(bian_temp1);
            if n<4
                app.bianLamp.Color="red";
                app.biantipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿHelp->ÿÿÿÿ";
                return;
            end
            if isempty(bian_temp1{1,2})
                app.bianLamp.Color="red";
                app.biantipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿHelp->ÿÿÿÿ";
                return;
            end
            bian_temp2=str2double(bian_temp1{1,4});
            app.bian=[];
            app.bian=[bian_temp1{1,1},bian_temp1{1,3},bian_temp2];
            index=isnan(bian_temp2);
            if length(index)<2
                app.bianLamp.Color="red";
                app.biantipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿHelp->ÿÿÿÿ";
                return;
            end
            app.A1EditField.Value=app.bian(index,1);
            app.A2EditField.Value=app.bian(index,2);
            app.bian(index,:)=[];
            app.bian=[app.bian,app.bian(:,3)];
            app.bianLamp.Color="green";
            app.biantipEditField.Value="ÿÿÿÿÿÿÿÿ";
            
            
            
            
        end

        % Button pushed function: jiaoButton
        function jiaoButtonPushed(app, event)
            if app.jiaosrcTextArea.Value==""
                app.jiaoLamp.Color="red";
                app.jiaotipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ";
                return;
            end
            
            value=strings;
            for i=1:length(app.jiaosrcTextArea.Value)
                value=value+app.jiaosrcTextArea.Value{i,1};
            end
            jiao_temp1 = textscan(value,'%c %f %f %f %f %c %f %c %f %c');
            
            [~,n]=size(jiao_temp1);
            if n<10
                app.jiaoLamp.Color="red";
                app.jiaotipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿHelp->ÿÿÿÿ";
                return;
            end
            if isempty(jiao_temp1{1,2})
                app.jiaoLamp.Color="red";
                app.jiaotipEditField.Value="ÿÿÿÿÿÿÿÿÿÿÿÿÿÿHelp->ÿÿÿÿ";
                return;
            end
            jiao_temp2=(jiao_temp1{:,5}+jiao_temp1{:,7}./60+jiao_temp1{:,9}./3600)/180*pi;
            app.jiao=[];
            app.jiao=[jiao_temp1{1,2},jiao_temp1{1,3},jiao_temp1{1,4},jiao_temp2,jiao_temp2];
            app.jiaoLamp.Color="green";
            app.jiaotipEditField.Value="ÿÿÿÿÿÿÿÿ";
            
        end

        % Value changed function: jiaosrcTextArea
        function jiaosrcTextAreaValueChanged(app, event)
            app.jiaosrcTextArea.Value = app.jiaosrcTextArea.Value;
            
        end

        % Menu selected function: Menu_8
        function Menu_8Selected(app, event)
            A = app.biansrcTextArea.Value;
            [file,path] = uiputfile('*.xlsx');
            if path==0
                return;
            end
            filename = fullfile(path,file);
            xlswrite(filename,A);
            app.biantipEditField.Value="ÿÿÿÿÿÿÿÿ"+filename;
        end

        % Menu selected function: Menu_9
        function Menu_9Selected(app, event)
            A = app.jiaosrcTextArea.Value;
            [file,path] = uiputfile('*.xlsx');
            if path==0
                return;
            end
            filename = fullfile(path,file);
            xlswrite(filename,A);
            app.jiaotipEditField.Value="ÿÿÿÿÿÿÿÿ"+filename;
        end

        % Value changed function: Spinner
        function SpinnerValueChanged(app, event)
             
            value = app.Spinner.Value;
            if value==0
                return;
            end
            if value>=app.count
                app.Spinner.Value=app.count;
                return;
            end
            app.UITableB.Data=app.data{value,1};
            app.UITableP.Data=app.data{value,2};
            app.UITablel.Data=app.data{value,3};
            app.UITable.Data=app.data{value,4};
            app.UITablex.Data=app.data{value,5};
            %app.UITableP.Data=P;
            %app.UITablel.Data=l;
        end

        % Menu selected function: BMenu
        function BMenuSelected(app, event)
            A = app.UITableB.Data;
            [file,path] = uiputfile('*.xlsx');
            if path==0
                return;
            end
            filename = fullfile(path,file);
            xlswrite(filename,A);
            app.TextArea.Value="ÿÿÿÿÿÿÿ"+filename;
        end

        % Menu selected function: PMenu
        function PMenuSelected(app, event)
            A = app.UITableP.Data;
            [file,path] = uiputfile('*.xlsx');
            if path==0
                return;
            end
            filename = fullfile(path,file);
            xlswrite(filename,A);
            app.TextArea.Value="ÿÿÿÿÿÿÿ"+filename;
        end

        % Menu selected function: lMenu
        function lMenuSelected(app, event)
            A = app.UITablel.Data;
            [file,path] = uiputfile('*.xlsx');
            if path==0
                return;
            end
            filename = fullfile(path,file);
            xlswrite(filename,A);
            app.TextArea.Value="ÿÿÿÿÿÿÿ"+filename;
        end

        % Menu selected function: VMenu
        function VMenuSelected(app, event)
            A = app.UITable.Data;
            [file,path] = uiputfile('*.xlsx');
            if path==0
                return;
            end
            filename = fullfile(path,file);
            xlswrite(filename,A);
            app.TextArea.Value="ÿÿÿÿÿÿÿ"+filename;
        end

        % Menu selected function: xMenu
        function xMenuSelected(app, event)
            A = app.UITablex.Data;
            [file,path] = uiputfile('*.xlsx');
            if path==0
                return;
            end
            filename = fullfile(path,file);
            xlswrite(filename,A);
            app.TextArea.Value="ÿÿÿÿÿÿÿ"+filename;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ErrorAjustment and hide until all components are created
            app.ErrorAjustment = uifigure('Visible', 'off');
            app.ErrorAjustment.Position = [100 100 1111 721];
            app.ErrorAjustment.Name = 'Measurement Ajustment for Traverse Network';
            app.ErrorAjustment.Scrollable = 'on';

            % Create HomeMenu
            app.HomeMenu = uimenu(app.ErrorAjustment);
            app.HomeMenu.Text = 'Home';

            % Create MenuSave
            app.MenuSave = uimenu(app.HomeMenu);
            app.MenuSave.MenuSelectedFcn = createCallbackFcn(app, @MenuSaveSelected, true);
            app.MenuSave.Text = 'ÿÿÿÿ';

            % Create Menu_2
            app.Menu_2 = uimenu(app.HomeMenu);
            app.Menu_2.MenuSelectedFcn = createCallbackFcn(app, @Menu_2Selected, true);
            app.Menu_2.Text = 'ÿÿÿÿ';

            % Create Menu_8
            app.Menu_8 = uimenu(app.HomeMenu);
            app.Menu_8.MenuSelectedFcn = createCallbackFcn(app, @Menu_8Selected, true);
            app.Menu_8.Text = 'ÿÿÿÿÿ';

            % Create Menu_9
            app.Menu_9 = uimenu(app.HomeMenu);
            app.Menu_9.MenuSelectedFcn = createCallbackFcn(app, @Menu_9Selected, true);
            app.Menu_9.Text = 'ÿÿÿÿÿ';

            % Create Menu_10
            app.Menu_10 = uimenu(app.HomeMenu);
            app.Menu_10.Text = 'ÿÿÿÿÿÿ';

            % Create BMenu
            app.BMenu = uimenu(app.Menu_10);
            app.BMenu.MenuSelectedFcn = createCallbackFcn(app, @BMenuSelected, true);
            app.BMenu.Text = 'ÿÿB';

            % Create PMenu
            app.PMenu = uimenu(app.Menu_10);
            app.PMenu.MenuSelectedFcn = createCallbackFcn(app, @PMenuSelected, true);
            app.PMenu.Text = 'ÿÿP';

            % Create lMenu
            app.lMenu = uimenu(app.Menu_10);
            app.lMenu.MenuSelectedFcn = createCallbackFcn(app, @lMenuSelected, true);
            app.lMenu.Text = 'ÿÿl';

            % Create VMenu
            app.VMenu = uimenu(app.Menu_10);
            app.VMenu.MenuSelectedFcn = createCallbackFcn(app, @VMenuSelected, true);
            app.VMenu.Text = 'ÿÿV';

            % Create xMenu
            app.xMenu = uimenu(app.Menu_10);
            app.xMenu.MenuSelectedFcn = createCallbackFcn(app, @xMenuSelected, true);
            app.xMenu.Text = 'ÿÿx';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.ErrorAjustment);
            app.HelpMenu.Text = 'Help';

            % Create Menu_4
            app.Menu_4 = uimenu(app.HelpMenu);
            app.Menu_4.MenuSelectedFcn = createCallbackFcn(app, @Menu_4Selected, true);
            app.Menu_4.Text = 'ÿÿÿÿ';

            % Create Menu_5
            app.Menu_5 = uimenu(app.HelpMenu);
            app.Menu_5.MenuSelectedFcn = createCallbackFcn(app, @Menu_5Selected, true);
            app.Menu_5.Text = 'ÿÿÿÿ';

            % Create AboutMenu
            app.AboutMenu = uimenu(app.ErrorAjustment);
            app.AboutMenu.Text = 'About';

            % Create Menu_6
            app.Menu_6 = uimenu(app.AboutMenu);
            app.Menu_6.MenuSelectedFcn = createCallbackFcn(app, @Menu_6Selected, true);
            app.Menu_6.Text = 'ÿÿÿÿ';

            % Create Menu_7
            app.Menu_7 = uimenu(app.AboutMenu);
            app.Menu_7.MenuSelectedFcn = createCallbackFcn(app, @Menu_7Selected, true);
            app.Menu_7.Text = 'ÿÿÿÿ';

            % Create computeButton
            app.computeButton = uibutton(app.ErrorAjustment, 'push');
            app.computeButton.ButtonPushedFcn = createCallbackFcn(app, @computeButtonPushed, true);
            app.computeButton.FontSize = 18;
            app.computeButton.FontWeight = 'bold';
            app.computeButton.FontColor = [0 0 1];
            app.computeButton.Position = [29 349 492 36];
            app.computeButton.Text = 'ÿÿÿÿ';

            % Create Label
            app.Label = uilabel(app.ErrorAjustment);
            app.Label.HorizontalAlignment = 'center';
            app.Label.FontSize = 16;
            app.Label.FontWeight = 'bold';
            app.Label.Position = [29 314 494 22];
            app.Label.Text = 'ÿÿÿÿ';

            % Create TextArea
            app.TextArea = uitextarea(app.ErrorAjustment);
            app.TextArea.FontSize = 14;
            app.TextArea.Position = [29 1 494 302];

            % Create EditField_3Label
            app.EditField_3Label = uilabel(app.ErrorAjustment);
            app.EditField_3Label.HorizontalAlignment = 'right';
            app.EditField_3Label.FontSize = 14;
            app.EditField_3Label.FontWeight = 'bold';
            app.EditField_3Label.Position = [27 431 75 22];
            app.EditField_3Label.Text = 'ÿÿÿÿÿ';

            % Create A1EditField
            app.A1EditField = uieditfield(app.ErrorAjustment, 'numeric');
            app.A1EditField.Editable = 'off';
            app.A1EditField.Position = [110 431 100 22];

            % Create EditField_4Label
            app.EditField_4Label = uilabel(app.ErrorAjustment);
            app.EditField_4Label.HorizontalAlignment = 'right';
            app.EditField_4Label.FontSize = 14;
            app.EditField_4Label.FontWeight = 'bold';
            app.EditField_4Label.Position = [28 393 75 22];
            app.EditField_4Label.Text = 'ÿÿÿÿÿ';

            % Create A2EditField
            app.A2EditField = uieditfield(app.ErrorAjustment, 'numeric');
            app.A2EditField.Editable = 'off';
            app.A2EditField.Position = [110 393 100 22];

            % Create XmEditFieldLabel
            app.XmEditFieldLabel = uilabel(app.ErrorAjustment);
            app.XmEditFieldLabel.HorizontalAlignment = 'right';
            app.XmEditFieldLabel.FontSize = 14;
            app.XmEditFieldLabel.FontWeight = 'bold';
            app.XmEditFieldLabel.Position = [232 431 37 22];
            app.XmEditFieldLabel.Text = 'X(m)';

            % Create A1xEditField
            app.A1xEditField = uieditfield(app.ErrorAjustment, 'numeric');
            app.A1xEditField.ValueDisplayFormat = '%.3f';
            app.A1xEditField.Position = [272 431 100 22];

            % Create YmEditFieldLabel
            app.YmEditFieldLabel = uilabel(app.ErrorAjustment);
            app.YmEditFieldLabel.HorizontalAlignment = 'right';
            app.YmEditFieldLabel.FontSize = 14;
            app.YmEditFieldLabel.FontWeight = 'bold';
            app.YmEditFieldLabel.Position = [381 431 37 22];
            app.YmEditFieldLabel.Text = 'Y(m)';

            % Create A1yEditField
            app.A1yEditField = uieditfield(app.ErrorAjustment, 'numeric');
            app.A1yEditField.ValueDisplayFormat = '%.3f';
            app.A1yEditField.Position = [421 431 100 22];

            % Create YmEditField_2Label
            app.YmEditField_2Label = uilabel(app.ErrorAjustment);
            app.YmEditField_2Label.HorizontalAlignment = 'right';
            app.YmEditField_2Label.FontSize = 14;
            app.YmEditField_2Label.FontWeight = 'bold';
            app.YmEditField_2Label.Position = [382 393 37 22];
            app.YmEditField_2Label.Text = 'Y(m)';

            % Create A2yEditField
            app.A2yEditField = uieditfield(app.ErrorAjustment, 'numeric');
            app.A2yEditField.ValueDisplayFormat = '%.3f';
            app.A2yEditField.Position = [421 393 100 22];

            % Create XmEditField_2Label
            app.XmEditField_2Label = uilabel(app.ErrorAjustment);
            app.XmEditField_2Label.HorizontalAlignment = 'right';
            app.XmEditField_2Label.FontSize = 14;
            app.XmEditField_2Label.FontWeight = 'bold';
            app.XmEditField_2Label.Position = [233 393 37 22];
            app.XmEditField_2Label.Text = 'X(m)';

            % Create A2xEditField
            app.A2xEditField = uieditfield(app.ErrorAjustment, 'numeric');
            app.A2xEditField.ValueDisplayFormat = '%.3f';
            app.A2xEditField.Position = [272 393 100 22];

            % Create UIAxes
            app.UIAxes = uiaxes(app.ErrorAjustment);
            title(app.UIAxes, 'ÿÿÿÿÿÿ ÿÿm')
            xlabel(app.UIAxes, 'Y')
            ylabel(app.UIAxes, 'X')
            app.UIAxes.TitleFontWeight = 'bold';
            app.UIAxes.Position = [599 335 456 381];

            % Create bianButton
            app.bianButton = uibutton(app.ErrorAjustment, 'push');
            app.bianButton.ButtonPushedFcn = createCallbackFcn(app, @bianButtonPushed, true);
            app.bianButton.FontWeight = 'bold';
            app.bianButton.Position = [48 494 199 24];
            app.bianButton.Text = 'ÿÿÿÿÿ';

            % Create jiaoButton
            app.jiaoButton = uibutton(app.ErrorAjustment, 'push');
            app.jiaoButton.ButtonPushedFcn = createCallbackFcn(app, @jiaoButtonPushed, true);
            app.jiaoButton.FontWeight = 'bold';
            app.jiaoButton.Position = [320 494 201 24];
            app.jiaoButton.Text = 'ÿÿÿÿÿ';

            % Create bianLamp
            app.bianLamp = uilamp(app.ErrorAjustment);
            app.bianLamp.Position = [28 496 20 20];
            app.bianLamp.Color = [1 0 0];

            % Create jiaoLamp
            app.jiaoLamp = uilamp(app.ErrorAjustment);
            app.jiaoLamp.Position = [301 496 20 20];
            app.jiaoLamp.Color = [1 0 0];

            % Create biantipEditField
            app.biantipEditField = uieditfield(app.ErrorAjustment, 'text');
            app.biantipEditField.Editable = 'off';
            app.biantipEditField.Position = [29 461 218 22];

            % Create jiaotipEditField
            app.jiaotipEditField = uieditfield(app.ErrorAjustment, 'text');
            app.jiaotipEditField.Editable = 'off';
            app.jiaotipEditField.Position = [306 461 215 22];

            % Create Label_2
            app.Label_2 = uilabel(app.ErrorAjustment);
            app.Label_2.HorizontalAlignment = 'center';
            app.Label_2.FontSize = 14;
            app.Label_2.FontWeight = 'bold';
            app.Label_2.Position = [301 694 220 22];
            app.Label_2.Text = 'ÿÿÿÿÿÿÿÿÿÿ';

            % Create jiaosrcTextArea
            app.jiaosrcTextArea = uitextarea(app.ErrorAjustment);
            app.jiaosrcTextArea.ValueChangedFcn = createCallbackFcn(app, @jiaosrcTextAreaValueChanged, true);
            app.jiaosrcTextArea.Position = [301 528 220 166];

            % Create Label_3
            app.Label_3 = uilabel(app.ErrorAjustment);
            app.Label_3.HorizontalAlignment = 'center';
            app.Label_3.FontSize = 14;
            app.Label_3.FontWeight = 'bold';
            app.Label_3.Position = [27 694 220 22];
            app.Label_3.Text = 'ÿÿÿÿÿÿÿÿÿÿ';

            % Create biansrcTextArea
            app.biansrcTextArea = uitextarea(app.ErrorAjustment);
            app.biansrcTextArea.Position = [27 528 220 166];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.ErrorAjustment);
            app.TabGroup.Position = [554 1 543 302];

            % Create BTab
            app.BTab = uitab(app.TabGroup);
            app.BTab.Title = 'B';

            % Create UITableB
            app.UITableB = uitable(app.BTab);
            app.UITableB.ColumnName = '';
            app.UITableB.RowName = {};
            app.UITableB.ColumnEditable = true;
            app.UITableB.Position = [1 1 541 276];

            % Create PTab
            app.PTab = uitab(app.TabGroup);
            app.PTab.Title = 'P';

            % Create UITableP
            app.UITableP = uitable(app.PTab);
            app.UITableP.ColumnName = '';
            app.UITableP.RowName = {};
            app.UITableP.ColumnEditable = true;
            app.UITableP.Position = [1 1 541 276];

            % Create lTab
            app.lTab = uitab(app.TabGroup);
            app.lTab.Title = 'l';

            % Create UITablel
            app.UITablel = uitable(app.lTab);
            app.UITablel.ColumnName = '';
            app.UITablel.RowName = {};
            app.UITablel.ColumnEditable = true;
            app.UITablel.Position = [2 1 540 276];

            % Create VTab
            app.VTab = uitab(app.TabGroup);
            app.VTab.Title = 'V';

            % Create UITable
            app.UITable = uitable(app.VTab);
            app.UITable.ColumnName = '';
            app.UITable.RowName = {};
            app.UITable.ColumnEditable = true;
            app.UITable.Position = [1 1 541 276];

            % Create xTab
            app.xTab = uitab(app.TabGroup);
            app.xTab.Title = 'x';

            % Create UITablex
            app.UITablex = uitable(app.xTab);
            app.UITablex.ColumnName = '';
            app.UITablex.RowName = {};
            app.UITablex.ColumnEditable = true;
            app.UITablex.Position = [1 1 541 276];

            % Create Label_4
            app.Label_4 = uilabel(app.ErrorAjustment);
            app.Label_4.HorizontalAlignment = 'right';
            app.Label_4.FontWeight = 'bold';
            app.Label_4.Position = [797 314 137 22];
            app.Label_4.Text = 'ÿÿÿÿÿÿÿÿÿÿÿ';

            % Create Spinner
            app.Spinner = uispinner(app.ErrorAjustment);
            app.Spinner.Limits = [0 1];
            app.Spinner.ValueDisplayFormat = '%.0f';
            app.Spinner.ValueChangedFcn = createCallbackFcn(app, @SpinnerValueChanged, true);
            app.Spinner.Position = [949 314 124 22];

            % Create Label_5
            app.Label_5 = uilabel(app.ErrorAjustment);
            app.Label_5.HorizontalAlignment = 'right';
            app.Label_5.FontWeight = 'bold';
            app.Label_5.Position = [581 314 113 22];
            app.Label_5.Text = 'ÿÿÿÿÿÿÿÿÿ';

            % Create EditFieldCount
            app.EditFieldCount = uieditfield(app.ErrorAjustment, 'numeric');
            app.EditFieldCount.ValueDisplayFormat = '%.0f';
            app.EditFieldCount.Editable = 'off';
            app.EditFieldCount.Position = [709 314 77 22];

            % Show the figure after all components are created
            app.ErrorAjustment.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = daoxianwang_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ErrorAjustment)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ErrorAjustment)
        end
    end
end