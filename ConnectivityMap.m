% Copyright (C) 2019 JT McCoy and L Auret

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.

% Controller maintenance prioritisation with LoopRank for a milling and flotation plant
% submitted to IFAC MMM 2019 Conference

% JT McCoy and L Auret
% Department of Process Engineering
% University of Stellenbosch

% Acknowledgements to Brian Lindner for supplying this code.
% For open source version, see: 
% https://github.com/ProcessMonitoringStellenboschUniversity/IFAC-Causality-Analysis

classdef ConnectivityMap<handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ConnGraph
        ExtractionMethod
        Symptoms
        FaultGraph
        Roots
        PathEdges
        PathNodes
    end
    
    methods
        
        function this  = ConnectivityMap()
            %CONNECTIVITYMAP Class constructor
            
        end
        
        function drawMap(this,cm,names,extrMeth)
           
%             f= figure;
            graph = this.constructGraph(cm,names,false);
            graph.view
            
%             set(0, 'ShowHiddenHandles', 'on')
%             bgf = gcf;
%             c1 = get(bgf, 'Children');
%             h1 = copyobj(c1(9), f);
%             close(bgf)
            % Package results and info
            this.ConnGraph = graph;
            this.ExtractionMethod = extrMeth;
        end
        
        function  comparison(this,cm1,cm2,refGraph)
            %COMPARISON compares two connectivity graphs
            
            M = length(cm1);
            names = get(refGraph.Nodes,'ID');
            sigInd1 = cm1>0;
            sigInd2 = cm2>0;
            
            compCM = zeros(M,M);
            compCM(sigInd1==1 & sigInd2==1) = 1; %no change
            change = sigInd2-sigInd1;
            compCM(change>0) = 2; %move above sig threshold
            compCM(change<0) = 3; %move below sig threshold
            
            compGraph = this.constructGraph(compCM,names,false);
            
            % Colour increased edges blue
            [r,c] = find(compCM==2);
            if ~isempty(r)
                sinkSource(:,1) = names(r);
                sinkSource(:,2) = names(c);
                for ss = 1:length(sinkSource)
                    incEdges(ss) = compGraph.getedgesbynodeid(sinkSource(ss,1),sinkSource(ss,2));
                end
                set(incEdges,'LineColor',[0.1 0.1 0.9])
            end
            
            % Colour mass balance edges Blue
            [r,c] = find(compCM==3);
            if ~isempty(r)
                sinkSource(:,1) = names(r);
                sinkSource(:,2) = names(c);
                for ss = 1:length(sinkSource)
                    decEdges(ss) = compGraph.getedgesbynodeid(sinkSource(ss,1),sinkSource(ss,2));
                end
                set(decEdges,'LineColor',[0.9 0.1 0.1])
            end
            
        end
        
        function connectivityUI(this,cm,names)
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here
            
            % Create a figure with size relative to screen size:
            % Using 2/3 of screen height, with 1/2 of screen width
            scrSize = get(groot,'ScreenSize');
            f = figure('Visible','on','Position',[scrSize(3)/2 scrSize(4)/6 scrSize(3)/2 scrSize(4)*2/3]);
           
            % Draw graph (find out how to do it with constructGraph
            % function)
            %                 cg = biograph(cm,names,'ShowArrows','on',...
            %                     'ShowWeights','off','NodeAutoSize','on');
            %                 set(cg.nodes,'FontSize',16);
            %                 set(cg.nodes,'LineWidth',1.5);
            %                 set(cg.nodes,'Color',[0.9 0.9 0.9]);
            %                 set(cg.nodes,'LineColor',[0 0 0])
            %                 set(cg.edges,'LineColor',[0 0 0])
            %                 set(cg.nodes,'Shape','circle')
            %                 dolayout(cg);
            cg = this.constructGraph(cm,names,0);
            cg.view
            set(0, 'ShowHiddenHandles', 'on')
            bgf = gcf;
            c1 = get(bgf, 'Children');
            h1 = copyobj(c1(9), f);
            
            close(bgf)
            m = 50;
            h1.Units = 'pixels';
            
            
            figWidth = f.Position(3)-150-m;
            scle = figWidth/h1.Position(3);
            figHeight = scle*h1.Position(4);
            if figHeight>f.Position(4)-2*m
                scle = (f.Position(4)-2*m)/figHeight;
                figHeight = f.Position(4)-2*m;
                figWidth = scle*figWidth;
            end
            h1.Position = [150,m,figWidth,figHeight];
            f.Units = 'normalized';
            
            % panel = uipanel(f);
            %             bg = uibuttongroup('Title','Options','Visible', 'off','Position',[0 0 .2 1]);
            hedit  = uicontrol('Style','edit','String','0','Tag','Edit',...
                'Position',[10 500 100 30],'Callback',{@thresh_Callback,cm,names,f});
            hslide  = uicontrol('Style','slider','Min',0,...
                'Max',max(max(cm)),'Value',0,'Tag', 'Slider',...
                'Position',[10 450 100 30],'Callback',{@thresh_Callback,cm,names,f});
            hstext = uicontrol('Style','text','String','Threshold = ',...
                'Position',[10 530 100 30],'Tag','SlideText');
            hcheck = uicontrol('Style','checkbox','String','Fix Positions',...
                'Position',[10 350 100 30],'Callback',{@checkbox_Callback,cg},...
                'Tag','check1');
            htext  = uicontrol('Style','text','String','Select Extraction Method',...
                'Position',[10 300 100 30]);
            hpopup = uicontrol('Style','popupmenu',...
                'String',{'Granger Causality','Transfer Entropy'},...
                'Position',[10 250 100 30],...
                'Callback',@popup_menu_Callback);
            
            hcalc   = uicontrol('Style','pushbutton','String','Calculate',...
                'Position',[10 150 100 30],...
                'Callback',{@calcbutton_Callback,cm,f});
            
            
            align([hedit,hslide,hstext,hcheck,htext,hpopup,hcalc],'Center','None');
            
            
            % Initialize the UI.
            % Change units to normalized so components resize automatically.
            hedit.Units = 'normalized';
            hslide.Units = 'normalized';
            hstext.Units = 'normalized';
            hcheck.Units = 'normalized';
            hcalc.Units = 'normalized';
            htext.Units = 'normalized';
            hpopup.Units = 'normalized';
            % Assign a name to appear in the window title.
            f.Name = 'Causality GUI';
            
            % Move the window to the center of the screen.
            % movegui(f,'center')
            
            % Make the UI visible.
            %             bg.Visible = 'on';
            f.Visible = 'on';
            
            
            %  Pop-up menu callback. Read the pop-up menu Value property to
            %  determine which item is currently displayed and make it the
            %  current data. This callback automatically has access to
            %  current_data because this function is nested at a lower level.
            function popup_menu_Callback(source,eventdata)
                % Determine the selected data set.
                str = source.String;
                val = source.Value;
                % Set current data to the selected data set.
                switch str{val};
                    case 'Transfer Entropy' % User selects Peaks.
                        
                    case 'Granger Causality' % User selects Membrane.
                        
                end
            end
            
            % Push button callbacks. Each callback plots current_data in the
            % specified plot type.
            
            function checkbox_Callback(source,eventdata,cg)
                
                val = source.Value;
                source.UserData = get(cg.Nodes,'Position');
                
            end
            
            function calcbutton_Callback(source,eventdata,cg,f)
                cg.view
                set(0, 'ShowHiddenHandles', 'on')
                bgfig = gcf;
                c = get(bgfig, 'Children');
                ha = copyobj(c(9), f);
                close(bgfig)
                mrgin = 50;
                % ha = axes('Units','pixels','Position',[f.Position(3)/2,mrgin,f.Position(3)/2-mrgin,f.Position(4)/2-mrgin]);
                ha.Units = 'pixels';
                ha.Position = [f.Position(3)/2,mrgin,f.Position(3)/2-mrgin,f.Position(4)-mrgin];
                f.Units = 'normalized';
            end
            
            function thresh_Callback(source,eventdata,cm,names,f)
                
                %                 hst = findobj('Tag','SlideText');
                %                 hst.String = ['Threshold = ' num2str(source.Value)];
                f.Units = 'pixels';
                
                % Slider and edit box both change the threshold, so this is
                % to ensure that the value from each is read and that when
                % one changes so does the other
                if strcmp(source.Style,'slider')
                    threshVal = source.Value;
                    hed = findobj('Tag','Edit');
                    hed.String = num2str(threshVal);
                elseif strcmp(source.Style,'edit')
                    threshVal = str2num(source.String);
                    hsl = findobj('Tag','Slider');
                    % If threshold value is greater than maximum, just set
                    % the slider value to maximum (similarly when less than
                    % min)
                    if threshVal>max(max(cm))
                        hsl.Value = max(max(cm));
                    elseif threshVal<0
                        hsl.Value = 0;
                    else
                        hsl.Value = threshVal;
                    end
                end
                cmSig = cm;
                cmSig(cm<threshVal) = 0;
                % Draw graph (find out how to do it with constructGraph
                % function)
                graph = biograph(cmSig,names,'ShowArrows','on',...
                    'ShowWeights','off','NodeAutoSize','on');
                set(graph.nodes,'FontSize',16);
                set(graph.nodes,'LineWidth',1.5);
                set(graph.nodes,'Color',[0.9 0.9 0.9]);
                set(graph.nodes,'LineColor',[0 0 0])
                set(graph.edges,'LineColor',[0 0 0])
                set(graph.nodes,'Shape','circle')
                dolayout(graph);
                
                hc = findobj('Tag','check1');
                fixPos = hc.Value;
                if fixPos
                    pos = hc.UserData;
                    for nde = 1:length(graph.Nodes)
                        set(graph.Nodes(nde),'Position',pos{nde})
                    end
                    dolayout(graph,'Pathsonly',true);
                end
                
                graph.view
                set(0, 'ShowHiddenHandles', 'on')
                bgfig = gcf;
                c = get(bgfig, 'Children');
                f.Children(end).delete; % careful with this hard-code number!
                ha = copyobj(c(9), f);
                %                 ha.NextPlot = 'replacechildren';
                close(bgfig)
                mrgin = 50;
                % ha = axes('Units','pixels','Position',[f.Position(3)/2,mrgin,f.Position(3)/2-mrgin,f.Position(4)/2-mrgin]);
                ha.Units = 'pixels';
                fWidth = f.Position(3)-150-mrgin;
                scl = fWidth/ha.Position(3);
                fHeight = scl*ha.Position(4);
                if fHeight>f.Position(4)-2*mrgin
                    scl = (f.Position(4)-2*mrgin)/fHeight;
                    fHeight = f.Position(4)-2*mrgin;
                    fWidth = scl*fWidth;
                end
                ha.Position = [150,mrgin,fWidth,fHeight]; f.Units = 'normalized';
                
            end
            
        end
        
        function backProp(this,symptNodes)
            
            %get faultnodes from contributions/connchange
            faultGraph = this.ConnGraph.deepCopy;
            for m = 1:length(faultGraph.nodes)
                varNames{m} = faultGraph.nodes(m).ID;
            end
            %Reset colours
            set(faultGraph.nodes,'Color',[0.9 0.9 0.9])
            % Highlight symptom nodes in blue
            set(faultGraph.getnodesbyid(symptNodes),'Color',[0.3 0.3 0.9])
            
            numSympt=length(symptNodes);
            M = length(faultGraph.Nodes);
            
            %find paths for all fault nodes
            for sNode = 1:numSympt
                ancNodes{sNode} = getancestors(faultGraph.getnodesbyid(symptNodes(sNode)),M);
                tempInd = zeros(1,M);
                for an=1:length(ancNodes{sNode})
                    % finding the Node numbers for this symptom nodes ancestor nodes
                    tempInd(an,:)=strcmp(ancNodes{sNode}(an).ID,varNames);
                end
                AncNodesInd(sNode,:)=sum(tempInd,1);%index of ancestor nodes for each fault node
            end
            
            commonAncestors = varNames(sum(AncNodesInd,1)==numSympt);
            if isempty(commonAncestors)
                warning('No common ancestors found')
                path = cell.empty;
                pathEdges = cell.empty;
            else
                commAncInd = find(sum(AncNodesInd,1)==numSympt); % node index
                set(faultGraph.getnodesbyid(commonAncestors),'Color',[0.9 0.2 0.2]);
                
                pathEdges = [];
                for sNode = 1:numSympt
                    symptInd = find(strcmp(symptNodes{sNode},varNames));
                    for ca =1:length(commonAncestors)
                        %path from common ancestor to sympyom node
                        
                        [~,path{ca,sNode},~]=shortestpath(faultGraph,commAncInd(ca),symptInd);
                        
                        for s=2:length(path{ca,sNode})
                            sourcenode=faultGraph.nodes(path{ca,sNode}(s-1)).ID;
                            sink=faultGraph.nodes(path{ca,sNode}(s)).ID;
                            pathEdges = getedgesbynodeid(faultGraph,sourcenode,sink);
                            set(pathEdges,'LineColor',[0.9 0.1 0.1]);
                            set(pathEdges,'LineWidth',2);
                        end
                    end
                end
            end
            % Package results
            this.Symptoms = symptNodes;
            this.FaultGraph = faultGraph;
            this.Roots = commonAncestors;
            this.PathEdges = pathEdges;
            this.PathNodes = path;
        end% Back propagation
        
    end
    
    methods (Static)
        
        function graph = constructGraph(cm,names,weightSize)
            %CONSTRUCTGRAPH Constructs connectivity graph from connectivity matrix
            
            graph = biograph(cm,names,'ShowArrows','on','ShowWeights','off','NodeAutoSize','on');
            set(graph.nodes,'FontSize',16);
            set(graph.nodes,'LineWidth',1.5);
            set(graph.nodes,'Color',[0.9 0.9 0.9]);
            set(graph.nodes,'LineColor',[0 0 0])
            set(graph.edges,'LineColor',[0 0 0])
            set(graph.nodes,'Shape','circle')
            
            if weightSize
                % Edge thickness relates to weights
                %Scale weights;
                lwMax = 4;
                lwMin = 1;
                for eNum = 1:length(graph.Edges)
                    lw = (lwMax-lwMin)*(graph.edges(eNum).Weight)+lwMin;
                    set(graph.Edges(eNum),'LineWidth',lw)
                end
            end
            dolayout(graph);
            
        end
              
    end
end

