function GUI_depthPitch_LRAUV

global guivars s4 p11 p12 p13 p14 sc21 sc22 sc23 p41 br b11 b12 b2 b3 eleCangle dzC qual

% Create pop-up1 menu
%--------------------------------------------------------------------------
uicontrol('Style', 'popup',...
    'String', guivars.names,...
    'Units','normalized',...            % [left bottom width height]
    'Position', [0.65 0.2 0.1 0.03],...
    'Callback', @pop_up1_Callback_depthPitch_LRAUV);



% Create pop-up2 menu
%--------------------------------------------------------------------------
stringlist = {'Display data only','Add malfunction','Add condtion'};
uicontrol('Style', 'popup',...
    'String', stringlist,...
    'Units','normalized',...            % [left bottom width height]
    'Position', [0.75 0.2 0.1 0.03],...
    'Callback', @pop_up2_Callback_depthPitch_LRAUV);


% Toggle button: select data (brush)
%--------------------------------------------------------------------------
tb1 = uicontrol('Style', 'togglebutton', 'String', 'Select data',...
    'Units','normalized','Position', [0.65 0.14 0.1 0.05],...
    'Value',0,'Callback', @togglebutton_Callback_depthPitch_LRAUV);


% Push button: clear data
%--------------------------------------------------------------------------
pb1 = uicontrol('Style', 'pushbutton', 'String', 'Clear selection',...
    'Units','normalized','Position', [0.65 0.08 0.1 0.05],...
    'Callback', @pushbutton_Callback_depthPitch_LRAUV);

% turn off functionality if MATLAB ver is earlier than R2014b
if verLessThan('matlab','8.4')
    set([tb1,pb1],'Enable', 'off');
end

% Slider1: select filtering criteria threshold - Elevator angle
%--------------------------------------------------------------------------
% Add a text uicontrol to label the slider.
sld1 = uicontrol('Style','text',...
    'Units','normalized',...            % [left bottom width height]
    'Position', [0.755 0.1625 0.195 0.0275],...
    'HorizontalAlignment','left',...
    'String','Elevator angle (deg) threshold:');
% Create slider
sld1tx = uicontrol('Style', 'slider',...
    'Min',0,'Max',20,'Value',eleCangle,...
    'Units','normalized',...            % [left bottom width height]
    'Position', [0.755 0.14 0.195 0.0275],...
    'Callback', @slider1_Callback_depthPitch_LRAUV);
% Display val in window
Vedit1=uicontrol('Style','edit',...
    'Units','normalized',...            % [left bottom width height]
    'Position', [0.955 0.1475 0.035 0.025],... [0.8575 0.1705 0.035 0.025]
    'String',num2str(eleCangle,3),...
    'Callback', @Vedit1_Callback_depthPitch_LRAUV);

set([sld1,sld1tx,Vedit1],'Enable', 'off');


% Slider2: select filtering criteria threshold - Depth rate
%--------------------------------------------------------------------------
% Add a text uicontrol to label the slider.
sld2 = uicontrol('Style','text',...
    'Units','normalized',...            % [left bottom width height]
    'Position', [0.755 0.1025 0.195 0.0275],...
    'HorizontalAlignment','left',...
    'String','Depth rate (m/s) threshold:');
% Create slider
sld2tx = uicontrol('Style', 'slider',...
    'Min',0,'Max',1.6,'Value',dzC,...
    'Units','normalized',...            % [left bottom width height]
    'Position', [0.755 0.08 0.195 0.0275],...
    'Callback', @slider2_Callback_depthPitch_LRAUV);
% Display val in window
Vedit2=uicontrol('Style','edit',...
    'Units','normalized',...            % [left bottom width height]
    'Position', [0.955 0.0875 0.035 0.025],... [0.8575 0.1105 0.035 0.025]
    'String',num2str(dzC,3),...
    'Callback', @Vedit2_Callback_depthPitch_LRAUV);

set([sld2,sld2tx,Vedit2],'Enable', 'off');




% CALLBACK FUNCTIONS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% pop-up1 callback
%--------------------------------------------------------------------------
    function pop_up1_Callback_depthPitch_LRAUV(source,callbackdata)
        
        val = get(source,'Value');
        str = get(source,'String');
        ind = strcmpi(guivars.names,str{val});
        
        
            
            % set(p41,'XData', guivars.time);           % time vector
            set(p41,'YData', guivars.data(ind,:));    % chosen dataset
            ylabel(s4,guivars.names{ind},'fontweight','bold','fontsize',16);
            if strcmpi(str{val},'Ele angle of attack (deg)')
                set(s4,'YLim',[-15 15])
            else
                set(s4,'YLimMode','auto')
            end
            
            %{
            if ~strcmp(str{val},'Mission')
            else
                
                p41 = plot(guivars.mission.time , guivars.mission.z,'linewidth',2);
                set(gca,'YDir','reverse')
                legend(guivars.mission.namelist,'location','best','Interpreter','none')
                dynamicDateTicks(gca)
                set(gca,'fontsize',11,'fontWeight','bold')
                ylabel(s4,guivars.names{ind},'fontweight','bold','fontsize',16)
                legend(guivars.mission.namelist,'location','best','Interpreter','none')
                grid on;
            end
            %}
        drawnow
        
    end


% pop-up2 callback
%--------------------------------------------------------------------------
    function pop_up2_Callback_depthPitch_LRAUV(source,callbackdata)
   
        val = get(source,'Value');
        str = get(source,'String');
        
        if strcmpi('Add malfunction',str{val})
            
            set(sc22,'Visible','on');   set(sc23,'Visible','off');
            set(p13,'Visible','on');    set(p14,'Visible','off');
            set([sld1,sld1tx,Vedit1],'Enable', 'off');
            set([sld2,sld2tx,Vedit2],'Enable', 'off');
            
        elseif strcmpi('Add condtion',str{val})
            
            set([sc22, sc23],'Visible','on')
            set([p13, p14],'Visible','on')
            
            if ~verLessThan('matlab','8.4')
                set([sld1,sld1tx,Vedit1],'Enable', 'on');
                set([sld2,sld2tx,Vedit2],'Enable', 'on');
            end
        else
            set([sc22, sc23],'Visible','off')
            set([p13, p14],'Visible','off')
            set([sld1,sld1tx,Vedit1],'Enable', 'off');
            set([sld2,sld2tx,Vedit2],'Enable', 'off');
        end
    end



% select data toggle button
%--------------------------------------------------------------------------
    function togglebutton_Callback_depthPitch_LRAUV( hObject, eventdata, handles )
        % hObject    handle to togglebutton1 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        
        % Hint: get(hObject,'Value') returns toggle state of togglebutton1
        button_state = get(hObject,'Value');
        if button_state == get(hObject,'Max')
            
            % when button is pressed:
            set(br,'Enable','on')       % turn brush on
            
        elseif button_state == get(hObject,'Min')
            
            % when button is un-pressed:
            set(br,'Enable','off')      % turn brush off
            
            % get brushed data index
            b11=logical(get(p11, 'BrushData'));
            b12=logical(get(p12, 'BrushData'));
            b2 = logical(get(sc21, 'BrushData'));
            b3 = logical(get(p41, 'BrushData'));
            
            % prep brushed data index if empty
            if isempty(b11)
                b11 = false(size(get(p11, 'XData')));
                b12 = false(size(get(p12, 'XData')));
            end
            
            if isempty(b2)
                b2 = false(size(get(sc21, 'XData')));
            end
            
            if isempty(b3)
                b3 = false(size(get(p41, 'XData')));
            end
            
            % apply brush to other subplots:
            brInd  = b11 | b12 | b2 | b3;
            brInd(isnan(get(p12, 'YData'))) = false;
            set(p11, 'BrushData', uint8(brInd))
            set(p12, 'BrushData', uint8(brInd))
            set(p13, 'BrushData', uint8(brInd))
            set(sc21, 'BrushData', uint8(brInd))
            set(p41, 'BrushData', uint8(brInd))
            
        end
    end

% clear selection push button
%--------------------------------------------------------------------------
    function pushbutton_Callback_depthPitch_LRAUV(hObject, eventdata, handles)
        % hObject    handle to pushbutton1 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        
        set(p11,'BrushData',zeros(1,length(b11)));
        set(p12,'BrushData',zeros(1,length(b12)));
        set(sc21,'BrushData',zeros(1,length(b2)));
        set(sc22,'BrushData',zeros(1,length(b2)));
        set(sc23,'BrushData',zeros(1,length(b2)));
        set(p41,'BrushData',zeros(1,length(b3)));
        if any(get(p13,'BrushData'))
            set(p13,'BrushData',zeros(size(get(p13, 'YData'))))
        end
        if any(get(p14,'BrushData'))
            set(p14,'BrushData',zeros(size(get(p14, 'YData'))))
        end
    end


% txtbox Vedit1
%--------------------------------------------------------------------------
 function Vedit1_Callback_depthPitch_LRAUV(source,callbackdata)
     
     eleCangle = str2double(get(source,'String'));      % Elevator angle threshold
     set( sld1tx, 'Value', eleCangle );                 % Re-position slider
     
     qual = (guivars.data(7,:)<-eleCangle & guivars.data(1,:)>-dzC) |...
            (guivars.data(7,:)>eleCangle & guivars.data(1,:)<dzC) |...
            (guivars.data(7,:)<0 & guivars.data(1,:)<-0.4) |...
             guivars.data(4,:) < -35;
        
        set(sc23,'XData',guivars.slider.dzbi(qual),'YData',guivars.slider.pitchInd(qual))
        set(p14,'XData',guivars.slider.timeInd(qual),'YData',guivars.slider.zInd(qual))
 end


% slider1
%--------------------------------------------------------------------------
    function slider1_Callback_depthPitch_LRAUV(source,callbackdata)
        
        eleCangle = get(source,'Value');      % Elevator angle threshold
        set(Vedit1,'String',num2str(eleCangle,3));
        
        % ele   = guivars.data(7,:);
        % dzbi  = guivars.data(1,:);
        % pitch = guivars.data(4,:);
        
        qual = (guivars.data(7,:)<-eleCangle & guivars.data(1,:)>-dzC) |...
            (guivars.data(7,:)>eleCangle & guivars.data(1,:)<dzC) |...
            (guivars.data(7,:)<0 & guivars.data(1,:)<-0.4) |...
             guivars.data(4,:) < -35;
        
        set(sc23,'XData',guivars.slider.dzbi(qual),'YData',guivars.slider.pitchInd(qual))
        set(p14,'XData',guivars.slider.timeInd(qual),'YData',guivars.slider.zInd(qual))
        
        
    end


% txtbox Vedit1
%--------------------------------------------------------------------------
 function Vedit2_Callback_depthPitch_LRAUV(source,callbackdata)
     
     dzC = str2double(get(source,'String'));      % Elevator angle threshold
     set( sld2tx, 'Value', dzC );                 % Re-position slider
     
     qual = (guivars.data(7,:)<-eleCangle & guivars.data(1,:)>-dzC) |...
            (guivars.data(7,:)>eleCangle & guivars.data(1,:)<dzC) |...
            (guivars.data(7,:)<0 & guivars.data(1,:)<-0.4) |...
             guivars.data(4,:) < -35;
        
        set(sc23,'XData',guivars.slider.dzbi(qual),'YData',guivars.slider.pitchInd(qual))
        set(p14,'XData',guivars.slider.timeInd(qual),'YData',guivars.slider.zInd(qual))
 end

% slider2 
%--------------------------------------------------------------------------
    function slider2_Callback_depthPitch_LRAUV(source,callbackdata)
        
        dzC = get(source,'Value');      % Elevator angle threshold
        set(Vedit2,'String',num2str(dzC,3));
        
        % ele   = guivars.data(6,:);
        % dzbi  = guivars.data(2,:);
        qual = (guivars.data(7,:)<-eleCangle & guivars.data(1,:)>-dzC) |...
            (guivars.data(7,:)>eleCangle & guivars.data(1,:)<dzC) |...
            (guivars.data(7,:)<0 & guivars.data(1,:)<-0.4) |...
             guivars.data(4,:) < -35;
        
        set(sc23,'XData',guivars.slider.dzbi(qual),'YData',guivars.slider.pitchInd(qual))
        set(p14,'XData',guivars.slider.timeInd(qual),'YData',guivars.slider.zInd(qual))
        
        
    end
end