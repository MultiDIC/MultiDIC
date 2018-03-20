function [Answer,Canceled] = inputsdlg(Prompt, Title, Formats, DefAns, Options)
%INPUTSDLG Enhanced input dialog box supporting multiple data types
% ANSWER = INPUTSDLG(PROMPT) creates a modal dialog box that returns user
% input for multiple prompts in the cell array ANSWER. PROMPT is a 1-D
% cell array containing the PROMPT strings.
%
% Alternatively, PROMPT can have up to 4 columns. The first column
% sppecifies the prompt string. The second column to specify the struct
% field names to output ANSWER as a structure. The third column specifies
% units (i.e., post-fix labels to the right of controls) to display. The
% fourth column specifies the tooltip string. The tooltip string is ignored
% for text type.
%
% INPUTSDLG uses UIWAIT to suspend execution until the user responds.
%
% ANSWER = INPUTSDLG(PROMPT,NAME) specifies the title for the dialog.
%
% Note that INPUTSDLG(PROMPT) & INPUTSDLG(PROMPT,NAME) are similar to the
% standard INPUTDLG function, except for the dialog layout.
%
% ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS) can be used to specify the type
% of parameters to display with FORMATS matrix of structures. The dimension
% of FORMATS defines how PROMPT items are laid out in the dialog box. For
% example, if PROMPT has 6 elements and the size of FORMATS is 2x3 then,
% the items are shown in 2 rows, 3 columns format. The items in PROMPT
% correspond to a horizontal traversal of FORMATS.
%
% The fields in FORMATS structure are:
%
%   type     - Type of control ['check',{'edit'},'list','range','text',
%                               'color','table','button','none']
%   style    - UI control type used. One of:
%               {'checkbox'},                for 'check' type
%               {'edit'}                     for 'edit' type
%               {'listbox','popupmenu','radiobutton','togglebutton'}
%                                            for 'list' type
%               {'slider'}                   for 'range' type
%               {'text'}                     for 'text' type
%               {'edit'}                     for 'color' type
%               {'pushbutton'}               for 'button' and 'color' types
%               {'table'}                    for 'table' type
%   format   - Data format: ['text','date','float','integer','logical',
%                            'vector','file','dir']
%   limits   - [min max] (see below for details)
%   required -  'on'   - control must have an answer
%              {'off'} - control may return empty answer
%   items    - Type 'edit', Format 'file': File flter spec
%              Type 'list': Selection items (cell of strings)
%              Type 'table': Names of columns (cell of strings)
%              Type 'range': Slider step size spec [minor major]
%   size     - [width height] in pixels. Set to <=0 to auto-size.
%   enable   - Defines how to respond to mouse button clicks, including which
%              callback routines execute. One of:
%               {'on'}      - UI control is operational.
%                'inactive' ?UI control is not operational, but looks the
%                             same as when Enable is on.
%                'off'      ?UI uicontrol is not operational and its image
%                             is grayed out.
%   margin  -  A scalar or 2-element vector specifying the margin between control
%              and its labels in pixels.
%   labelloc - Prompt label location:
%               {'lefttop'}   - left of control, aligned to top
%                'leftmiddle' - left of control, aligned to middle
%                'leftbottom' - left of control, aligned to bottom
%                'topleft'    - above control, aligned left
%                'topcenter'  - above control, aligned center
%                'topright'   - above control, aligned right
%   unitsloc - Units label location:
%               {'righttop'}     - right of control, aligned to top
%                'rightmiddle'   - right of control, aligned to middle
%                'rightbottom'   - right of control, aligned to bottom
%                'bottomleft'    - below control, aligned left
%                'bottomcenter'  - below control, aligned center
%                'bottomright'   - below control, aligned right
%   callback - Defines callback funcion, a routine that executes whenever
%              you activate the uicontrol object. For the controls with
%              separate dialog, their callback functions are executed after
%              the dialog is closed. The callback function must be given as
%              a function handle with following syntax:
%
%                 my_callbackfcn(hobj,evt,handles,k)
%
%              where hobj and evt are the passed-through standard MATLAB
%              callback arguments, handles is a Nx3 array of dialog
%              objects. Here, the n-th row corresponds to the n-th PROMPT,
%              and handles(n,1) is the calling object handle (i.e., same as
%              hobj). handles(n,2) are the prompt texts and handles(n,3)
%              are the prompt unit texts.
%
%              For example, Formats(n,m).callback.ButtonDownFcn sets the
%              the button-press callback function.
%   span     - Defines size of objects in fields [rows columns]
%
% A missing field (either missing from FORMATS struct or the field value is
% left empty for an element) will be filled with a default field value.
%
% FORMATS type field defines what type of prompt item to be shown.
%
%   type  Description
%   -------------------------------------------------------------------
%   edit     Standard edit box (single or multi-line mode)
%   check    Check box for boolean item
%   list     Chose from a list of items ('listbox' style allows multiple item
%            selection)
%   range    Use slider to chose a value over a range
%   text     Static text (e.g., for instructions)
%   color    Color selection using 'uisetcolor'
%   button   Execute function defined in 'callback'
%   table    Uitable
%   none     A placeholder. May be used for its neighboring item to extend
%            over multiple columns or rows (i.e., "to merge cells")
%
% The allowed data format depends on the type of the field:
%
%   type    allowed format
%   --------------------------------------------
%   check     {logical}, integer, text
%   edit      {text}, date, float, integer, file, dir, vector
%   list      {integer}, text
%   range     {float}
%   color     {float}, integer
%   table     (a cell string to specify ColumnFormat of uiTable)
%   button    any data format allowed
%
% Formats 'file' and 'dir' for 'edit' type uses the standard UIGETFILE,
% UIPUTFILE, and UIGETDIR functions to retrieve a file or directory name.
%
% The role of limits field varies depending on other parameters:
%
%   style         role of limits
%   ---------------------------------------------------
%   checkbox      If data format is integer, limits(1) is the ANSWER value
%                 if the check box is not selected box is not selected and
%                 limits(2) is the ANSWER if the check box is selected. If
%                 data format is text, its limits field must be given as a
%                 cellstring array.
%   edit:text
%                 If diff(limits)>0, text aligns with the prompt label. If
%                 diff(limits)<0, tet aligns with the control.
%   edit::date
%                 limits must be a free-format date format string or a
%                 scalar value specifying the date format. Supported format
%                 numbers are: 0,1,2,6,13,14,15,16,23. The default date
%                 format is 2 ('mm/dd/yy'). See the tables in DATESTR help
%                 for the format definitions. As long as the user entry is
%                 a valid date/time expression, the dialog box
%                 automatically converts to the assigned format.
%   edit::float, edit::integer
%                 This style defines the range of allowed values if numeric
%                 vector is given, including the values specified. For
%                 other cases, use cell array to customize its
%                 validateattribute call. For example, {'positive','<',1}
%                 to specify value between 0 and 1, excluding 0 & 1.
%   edit::vector
%                 limits specifies the allowed number of elements in a
%                 vector. Use cell array to customize its validateattribute
%                 call if vector length is fixed or completely arbitrary.
%                 Note that 'column' attribute is always forced unless
%                 'row' attribute is explicitly specified.
%   edit::file
%                 If 0<diff(limits)<=1 uses UIGETFILE in single select
%                 mode with single-line edit. If diff(limits)>1 uses
%                 UIGETFILE in multi-select mode with multi-line edit. If
%                 diff(limits)<=0 usees UIPUTFILE with single-line edit
%   list::listbox If diff(limits)>1, multiple items can be selected. If
%                 auto-height and limits(1)>0, at most limits(1) lines will
%                 be shown.
%   slider        limits(1) defines the smallest value while
%                 limits(2) defines the largest value
%   color         'float' format limit must be [0 1] and 'integer' format
%                 limit must be [0 255]. If not, the behavior is not
%                 determined.
%   table         cell array defining ColumnWidths
%   none          If diff(limits)==0 space is left empty (default)
%                 If diff(limits)>0 : lets the item from left to extend
%                 If diff(limits)<0 : lets the item from above to extend
%                 NOTE: Use 'span' field of a control to automatically set
%                 the extension mode.
%
% Similar to how PROMPT strings are laid out, when FORMATS.style is set to
% either 'radiobutton' or 'togglebutton', FORMATS.items are laid out
% according to the dimension of FORMATS.items.
%
% There are two quick format options as well:
%
%  Quick Format Option 1 (mimicing INPUTDLG behavior):
%   FORMATS can specify the number of lines for each edit-type prompt in
%   FORMATS. FORMATS may be a constant value or a column vector having
%   one element per PROMPT that specifies how many lines per input field.
%   FORMATS may also be a matrix where the first column specifies how
%   many rows for the input field and the second column specifies how
%   many columns wide the input field should be.
%
%  Quick Format Option 2:
%   FORMATS can specify the types of controls and use their default
%   configurations. This option, however, cannot be used to specify
%   'list' control as its items are not specified. To use this option,
%   provide a string (if only 1 control) or a cell array of strings. If
%   a cell array is given, its dimension is used for the dialog
%   layout.
%
% ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS,DEFAULTANSWER) specifies the
% default answer to display for each PROMPT. For a non-tiled layout,
% DEFAULTANSWER must contain the same number of elements as PROMPT (that
% are not of 'none' style). If PROMPT does not provide ANSWER structure
% fields, DEFAULTANSWER should be a cell array with element type
% corresponding to FORMATS.format. Leave the cell element empty for a
% prompt with 'text' type. If ANSWER is a structure, DEFAULTANSWER must be
% a struct with the specified fields. (If additional fields are present in
% DEFAULTANSWER, they will be returned as parts of ANSWER.)
%
% For edit::file controls, a default answer that does not correspond to an
% existing file will be used as a default path and/or file name in the
% browse window.  It is passed as the DefaultName parameter to UIGETFILE or
% UIPUTFILE.
%
% To enable Tiled Mode, FORMATS must be given as a vector and DEFAULTANSWER
% must be given as a cell matrix or a struct vector. If FORMATS is a row
% vector, the dialog controls are tiled vertically; conversely if it is a
% column vector, the controls are tiled horizontally. If DEFAULTANSWER is
% given as a cell matrix, the number of rows of DEFAULTANSWER must match
% the number of elements of PROMPT. Each column of DEFAULTANSWER forms a
% tile row/column. If DEFAULTANSWER is given as a struct vector, each
% struct element forms a tile row/column.
%
% Empty rows and empty columns are automatically eliminated by default. To
% add an empty row or column, set the FORMATS entry of one of its cells to
% type = 'none' with its size = [H W] or [-Hmin -Wmin] specified to the
% desired width and height.
%
% ANSWER = INPUTSDLG(PROMPT,NAME,FORMATS,DEFAULTANSWER,OPTIONS) specifies
% additional options. If OPTIONS is the string 'on', the dialog is made
% resizable. If OPTIONS is a structure, the fields recognized are:
%
%  Option Field Description {} indicates the default value
%  ----------------------------------------------------------------------
%  Resize        Make dialog resizable: 'on' | {'off'}
%  WindowStyle   Sets dialog window style: {'normal'} | 'modal'
%  Interpreter   Label text interpreter: 'latex' | {'tex'} | 'none'
%  CancelButton  Show Cancel button: {'on'} | 'off'
%  ApplyButton   Adds Apply button: 'on' | {'off'}
%  Sep           Space b/w prompts in pixels: {10}
%  ButtonNames   Customize OK|Cancel|Apply button names: {up to 3 elements}
%  AlignControls Align adjacent controls in the same column: 'on' | {'off'}
%  FontSize      Customize font size. Default: get(0,'DefaultUicontrolFontSize')
%  CreateFcn     Callback function executed right after dialog creation
%                with syntax my_createfcn(hobj,evt,handles) with
%                standard MATLAB callback arguments, hobj & evt, and Nx3
%                array of handles. Here, the n-th row corresponds to the
%                n-th PROMPT, and handles(n,1) is the calling object handle
%                (i.e., same as hobj). handles(n,2) are the prompt texts
%                and handles(n,3) are the prompt unit texts.
%  DeleteFcn     Callback function executed just before deleting dialog.
%                The function syntax is the same as CreateFcn.
%
% [ANSWER,CANCELED] = INPUTSDLG(...) returns CANCELED = TRUE if user
% pressed Cancel button, closed the dialog, or pressed ESC. In such event,
% the content of ANSWER is set to the default values.
%
% Note on Apply Button feature. Pressing the Apply button makes the current
% change permanent. That is, pressing Cancel button after pressing Apply
% button only reverts ANSWER back to the states when the Apply button was
% pressed last. Also, if user pressed Apply button, CANCELED flag will not
% be set even if user canceled out of the dialog box.
%
% Examples:
%
% prompt={'Enter the matrix size for x^2:';'Enter the colormap name:'};
% name='Input for Peaks function';
% formats(1) = struct('type','edit','format','integer','limits',[1 inf]);
% formats(2) = struct('type','edit','format','text','limits',[0 1]);
% defaultanswer={20,'hsv'};
%
% [answer,canceled] = inputsdlg(prompt,name,formats,defaultanswer);
%
% formats(2).size = -1; % auto-expand width and auto-set height
% options.Resize='on';
% options.WindowStyle='normal';
% options.Interpreter='tex';
%
% answer = inputsdlg(prompt,name,formats,defaultanswer,options);
%
% prompt(:,2) = {'Ndim';'Cmap'};
% defaultanswer = struct(defaultanswer,prompt(:,2),1);
%
% answer = inputsdlg(prompt,name,formats,defaultanswer,options);
%
% See also INPUTDLG, DIALOG, ERRORDLG, HELPDLG, LISTDLG, MSGBOX,
%  QUESTDLG, UIGETFILE, UIPUTFILE, UIGETDIR, DATESTR, VALIDATEATTRIBUTE.

% Version 2.2.0 (June 25, 2015)
% Written by: Takeshi Ikuma
% Contributors: Andreas Greuer, Luke Reisner, Florian Hatz
% Created: Nov. 16, 2009
% Revision History:
%  v.1.1 (Nov. 19, 2009)
%  * Fixed bugs (reported by AG):
%   - not returning Canceled output
%   - erroneous struct output behavior
%   - error if all row elements of a column are auto-expandable
%  * Added Apply button option
%  * Added support for Units (label to the right of controls)
%  * Updated the help text
%  v.1.11 (Nov. 20, 2009)
%  * Fixed bugs (reported by AG):
%   - incorrect Canceled output when Cancel button is pressed
%  v.1.12 (Nov. 20, 2009)
%  * Fixed bugs (reported by AG):
%   - again incorrect Canceled output behavior
%  v.1.2 (May 20, 2010)
%  * Fixed bugs (reported by AG & Jason):
%   - Apply button->Canel button does not revert back to post-apply answers.
%   - Line 265 handles.Figure -> handles.fig
%  * Added edit::date support
%  * Added formats.enable support
%  * Added options.CancelButton support
%  * Added options.ButtonNames support
%  v.1.2.1 (June 11, 2010)
%  * Fixed default option bug (reported by Jason)
%  v.1.2.2 (July 15, 2010)
%  * Rewritten checkoptions() (to correct issues reported by Jason)
%  * Bug Fix: file & dir control enable config were interpreted backwards
%  v.1.2.3 (July 19, 2010)
%  * checkoptions() bug fix (to correct issues reported by Kevin)
%  v.1.3 (August 13, 2010, by Luke Reisner)
%  * Improved dialog layout:
%   - Less wasted space, better control distribution, more consistent margins
%   - Buttons are right-aligned per OS standards
%  * Changed edit::date to return a simple date vector (see DATEVEC help)
%  * Added support for free-form date format specifiers to edit::date
%  * Added ability to limit the number of displayed lines for a listbox
%  * Added ability to set default browse path/filename for edit::file controls
%  * Added options.AlignControls to align adjacent controls in the same column
%  * Added options.UnitsMargin to control spacing between controls and units
%  * Fixed bugs:
%   - Flickering or misplaced controls when dialog first appears
%   - Radiobutton and togglebutton controls couldn't be disabled
%   - Edit::integer controls allowed non-integer values
%   - Slider controls didn't auto-size properly
%   - Other minor miscellaneous bugs
%  v.2.0 (July 17, 2013, by T. Ikuma & F. Hatz)
%  * PROMPT(:,4) to specify tooltip strings
%  * Enabled Tiled Mode with DEFAULTANSWER & FORMATS dimension specs. See
%    inputsdlg_demo_struct
%  * Added types: table, color, button
%  * Added formats: logical and vector
%  * Added 'text' format for 'list' type
%  * Added 'callback' format field
%  * Added 'CreateFcn' option field
%  * Added 'DeleteFcn' option field
%  * Added 'required' format field
%  * Added 'labelloc' and 'unitsloc' format fields to customize location of
%    the labels
%  * removed UnitMargin option field and added 'margin' format field.
%  * Added 'span' format field to make spanning across multiple rows and
%    columns simpler
%  * Improved handling of 'Formats' (new variable 'span' / see inputsdlg_demo)
%  * A dialog with non-editable text only displays 'OK' button and does not
%    return any argument
%  * edit:file: diff(formats.limits)==0 => uiputfile
%  v.2.0.1 (Aug 05, 2013, by T. Ikuma)
%  * Bug fix on inputdlg compatible calls
%  v.2.0.2 (Aug 06, 2013, by T. Ikuma)
%  * Version compatibility fix for pre-R2013a
%  v.2.0.3 (Sep 07, 2013)
%  * Bug fix on parsing popup menu callback (reported by E. Morales)
%  v.2.0.4 (Nov 25, 2013)
%  * Bug fix on parsing empty callback format field (reported by David v.B.)
%  v.2.0.5 (Mar 28, 2014)
%  * Improved error message on FORMATS-PROMPT size mismatch.
%  * Fixed several bugs relating to neglected ToolTipString support
%    (thanks David v.B!)
%  * FORMATS.items for Slider style now sets uicontrol's SliderStep
%    property
%  v.2.0.6 (Jul 01, 2014)
%  * bug fix in the text positioning (~Line 2350)
%  * bug fix in check_formats to insert type='edit' as default
%  v.2.1.0 (Jul 03, 2014)
%  * Added support for 'none' Format type with positive size to act as a
%    placeholder for a row or column.
%  * Removed restriction on Cancel & Apply buttons when there are only text
%    controls present.
%  * Bug fix on setting default answer for list:radiobutton and
%    list:togglebutton
%  * Bug fix on ANSWER reporting non-empty default answer for
%    list:radiobutton and list:togglebutton with required='off'.
%  * Bug fix on multi-select file formats
%  v.2.1.1 (Sept 17, 2014)
%  * Improved backward compatibility
%  * Bug fix in setting default value for list:integer:radiobutton control
%  v.2.1.2 (Sept 19, 2014)
%  * Improved check_formats routine (faster and better backward
%    compatibility)
%  v.2.2 (June 25, 2015)
%  * Added logical:text control with its limits format field specifying the
%    returning string values
%  * Numeric editboxes uses validateattribute to check its value. Its
%    formats.limits field may now be set to cell array to specify the
%    attribute
%  * Fixed errorneous behavior when pressed Enter key in Edit control
%    containing invalid value
%  v.2.3 (June 26, 2015)
%  * Added the row vector option for edit:vector via specifying 'row' in
%    the cell-array limits definition
%  * Bug fix in edit:vector limits parser
%  v.2.3.1 (June 29, 2015)
%  * Bug fixes in default value processing
%  v.2.3.2 (June 30, 2015)
%  * Bug fix for list:text default value processing (case insensitive)

% to-do list
% * Add edit:font for the built-in uisetfont dialog
% * Support for DefaultXXX option fields to set default Formats field values
% * Auto-layout given figure's desired width
% * Limits-based option for positioning of text-type
% * Better spanning support

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% # of argument Check %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
   narginchk(0,5);
   nargoutchk(0,2);
catch
   error(nargchk(0,5,nargin)); %#ok
   error(nargchk(0,2,nargout)); %#ok
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle Input Args %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1, Prompt={}; end
if nargin<2, Title = ''; end
if nargin<3, Formats=struct([]); end
if nargin<4, DefAns = {}; end
if nargin<5, Options = struct([]); end

% Check Prompt input
[Prompt,FieldNames,Units,TooltipStrings,err] = checkprompt(Prompt);
if ~isempty(err), error(err{:}); end
NumQuest = numel(Prompt); % number of prompts

if isempty(Title)
   Title = ' ';
elseif iscellstr(Title)
   Title = Title{1}; % take the first entry
elseif ~ischar(Title)
   error('inputsdlg:InvalidInput','Title must be a string of cell string.');
end

% make sure that the Options is valid
[Options,FormatDefaultFields,err] = checkoptions(Options);
if ~isempty(err), error(err{:}); end

% make sure that the Formats structure is valid & fill it in default values
% as needed
[Formats,err] = checkformats(Formats,NumQuest,FormatDefaultFields);
if ~isempty(err), error(err{:}); end

% make sure that the DefAns is valid & set Answer using DefAns and default
% values if DefAns not given
[DefaultAnswer,TileMode,err] = checkdefaults(DefAns,Formats,FieldNames);
if ~isempty(err), error(err{:}); end
Answer = DefaultAnswer;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Dialog GUI %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% lay contents out on a dialog box
[Formats,handles,sinfo] ...
   = buildgui(Title,Prompt,Units,TooltipStrings,FieldNames,Formats,TileMode,Options);

% fill in the default answers, setup callbacks,
initgui();

IsRequired = arrayfun(@(fmt)strcmp(fmt.required,'on'),Formats);
DefaultReqMet = checkReq(DefaultAnswer);
ReqMet = false; % also modified by the Apply button

ReturnPress = false;
Applied = false; % set true by pressing Apply Button
Canceled = ~ishghandle(handles.fig);

% Go into uiwait if the figure handle is still valid.
% This is mostly the case during regular use.
try
   while ~(Canceled || ReqMet) % exit only if Canceled or all the rerequired fields are filled
      
      % Wait till uiresume is called
      uiwait(handles.fig);
      
      figflag = get(handles.fig,'UserData');
   
      % if Return key press released from uiwait, make sure associated
      % callback did not cause an error
      if ReturnPress
         ReturnPress = false;
         if strcmp(figflag,'Error')
            ReqMet = false;
            set(handles.fig,'UserData','');
            continue;
         end
      end
      
      % Check handle validity again since figure could be deleted externally
      Canceled = strcmp(figflag,'Cancel');
      
      if Canceled % return the default answer
         Answer = DefaultAnswer; % revert back to the default answer
      else
         Answer = getAnswer(); % get the final answers
         ReqMet = checkReq(Answer);
         if ~ReqMet
            h = errordlg('All required parameters must be filled.','Missing Required Value(s)','modal');
            uiwait(h);
         end
      end
   end
   
   % if user deletefcn defined, call it now
   if ~isempty(Options.DeleteFcn)
      Options.DeleteFcn(handles.fig,[],handles.ctrls);
   end
   
   % Close the figure if it's still open
   delete(handles.fig);
   
catch ME
   if ishghandle(handles.fig)
      delete(handles.fig);
      ME.getReport
      throw(ME);
   else
      error('Inputsdlg dialog window was closed externally');
   end
end

% If Canceled, convert Canceled to integer depending on Applied condition
Canceled = Canceled * (Canceled + (Applied && DefaultReqMet));
% 0 - OK pressed
% 1 - Canceled
% 2 - Canceled but prior to it Applied button has been pressed, filled in
%     all the required answers

% If Tiled, reshape the answer
Answer = reshape(Answer,NumQuest,max(TileMode));

% If FieldNames given, convert Answer to struct
Answer = selectivecell2struct(Answer,FieldNames);

% If dialog contains only non-editable texts, return w/o output argument
if nargout==0 && all(sinfo.istext) 
   clear Answer
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function initgui
      
      fig = handles.fig;
      
      % set OK button as the default
      fh = handle(fig);
      fh.setDefaultButton(handles.btns(1));
      
      % Set callbacks of figure/control buttons
      set(fig, 'UserData', 'Cancel',...
         'WindowKeyPressFcn', @doFigureKeyPress,...
         'ResizeFcn', @(~,~)resizegui(handles,sinfo,Options),...
         'CloseRequestFcn', @(hObj,evd)doOKCancel(hObj,evd,false));
      for k = 1 : numel(handles.btns)
         hbtn = handles.btns(k);
         switch get(hbtn, 'UserData')
            case 'OK'
               set(hbtn, 'Callback', @(hObj,evd)doOKCancel(hObj,evd,true));
            case 'Cancel'
               set(hbtn, 'Callback', @(hObj,evd)doOKCancel(hObj,evd,false));
            case 'Apply'
               set(hbtn, 'Callback', @doApply);
         end
      end
      
      % Set callback functions and set default values
      for k = 1:numel(Answer)
         
         h = handles.ctrls(k,1);
         val = Answer{k};
         fmt = Formats(k);
         cbnames = fieldnames(fmt.callback);
         
         ena = {'Enable',fmt.enable};
         fcnname = 'Callback';
         fcn = {};
         aux = {};
         
         % get callback function handles for custom callbacks to be
         % called from inputsdlg's default callback
         [tf,I] = ismember('callback',lower(cbnames));
         if tf
            cbfcn = fmt.callback.(cbnames{I});
         else
            cbfcn = {};
         end
         [tf,I] = ismember('buttondownfcn',lower(cbnames));
         if tf
            bdfcn = fmt.callback.(cbnames{I});
         else
            bdfcn = {};
         end
         
         switch fmt.style
            case 'checkbox'
               ansname = 'Value';
               if strcmp(fmt.format,'text')
                  val = strcmp(val,fmt.limits{2});
               end
            case 'edit'
               ansname = 'String';
               switch fmt.format
                  case {'integer','float'}
                     % for numeric edit box, check for the range & set mouse down behavior
                     fcn = @(hObj,evd)checkNumericRange(hObj,evd,k,cbfcn);
                     %aux = {'UserData',val}; % save the numeric data
                     val = num2str(val);
                  case 'date'
                     fcn = @(hObj,evd)checkDate(hObj,evd,k,cbfcn);
                  case 'vector'
                     % for vector edit box, check for the range & set mouse down behavior
                     fcn = @(hObj,evd)checkVector(hObj,evd,k,cbfcn);
                     
                     val = num2str(val);
                  case 'file'
                     mode = diff(fmt.limits);
                     fcnname = 'ButtonDownFcn';
                     if strcmp(ena{2},'on')
                        ena{2} = 'inactive';
                        fcn = @(hObj,evd)openFilePrompt(hObj,evd,k,bdfcn);
                     end
                     
                     val = cellstr(val);
                     if ~isempty(val)
                        dirname = fileparts(val{1});
                        if ~isdir(dirname), dirname = ''; end
                     else
                        dirname = '';
                     end
                     aux = {'UserData',dirname};
                     
                     if mode <= 1 % single-file
                        val = val{1};
                     end
                  case 'dir'
                     fcnname = 'ButtonDownFcn';
                     if strcmp(ena{2},'on')
                        ena{2} = 'inactive';
                        fcn = @(hObj,evd)openDirPrompt(hObj,evd,k,bdfcn);
                     end
               end
            case 'pushbutton'
               if strcmp(fmt.type,'color')
                  if strcmp(ena{2},'on')
                     fcn = @(hObj,evd)openColorPrompt(hObj,evd,k,bdfcn);
                  end
                  ansname = 'BackgroundColor';
               else
                  ansname = 'UserData';
               end
            case {'radiobutton', 'togglebutton'}
               fcnname = 'SelectionChangeFcn';
               ansname = 'SelectedObject';
               hbtn = get(h,'UserData'); % hButtons
               if strcmp(fmt.format,'integer')
                  val = hbtn(val);
               else
                  val = findobj(hbtn,'flat','String',val);
               end
               set(hbtn,ena{:});
               ena = {};
            case 'table'
               ansname = 'Data';
               fcnname = 'CellEditCallback';
            case 'slider'
               if strcmp(fmt.format,'integer')
                  fcn = @(hObj,evd)forceInteger(hObj,evd,k,cbfcn);
               end
               ansname = 'Value';
            case {'listbox' 'popupmenu'}
               ansname = 'Value';
         end
         
         % Set control's properties
         if ~strcmp(fmt.type,'text')
            set(h,ena{:},fcnname,fcn,ansname,val,aux{:});
            
            % Set custom callbacks
            if ~isempty(fmt.callback)
               for n = 1:numel(cbnames)
                  % set if the specfied callback function not already assigned
                  % as the part of INPUTSDLG functionality.
                  if ~strcmpi(cbnames{n},fcnname) || isempty(fcn)
                     cbfcn = @(hobj,evt)fmt.callback.(cbnames{n})(hobj,evt,handles.ctrls,k);
                     try
                        set(h,cbnames{n},cbfcn);
                     catch
                        error('Invalid callback function name.');
                     end
                  end
               end
            end
         end
      end
      
      % make sure we are on screen
      movegui(handles.fig)
      
      % if there is a figure out there and it's modal, we need to be modal too
      if ~isempty(gcbf) && strcmp(get(gcbf,'WindowStyle'),'modal')
         set(handles.fig,'WindowStyle','modal');
      end
      
      % if user createfcn defined, call it now
      if ~isempty(Options.CreateFcn)
         Options.CreateFcn(handles.fig,[],handles.ctrls);
      end
      
      set(handles.fig,'Visible','on');
      drawnow;
      
      % set focus on the first uicontol
      h = findobj(handles.ctrls(:,1),'flat','-not','Style','text');
      if ~isempty(h)
         h = h(1);
         switch get(h,'type')
            case 'uicontrol', uicontrol(h);
            case {'uipanel' 'uibuttongroup'}
               hsel = get(h,'SelectedObject');
               if isempty(hsel), hsel = get(h,'Children'); hsel(1:end-1) = []; end
               uicontrol(hsel);
            case 'uitable', uitable(h);
         end
      end
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function answer = getAnswer()
      
      answer = cell(size(Answer));
      
      % retrieve answer from controls
      for i = 1:numel(answer)
         
         h = handles.ctrls(i,1);
         fmt = Formats(i);
         
         switch fmt.style
            case 'checkbox'
               val = get(h,'Value');
               
               if ~isempty(val)
                  switch fmt.format
                     case 'text' % 'popupmenu' 'listbox'
                        answer(i) = fmt.limits(val+1);
                     case 'logical' % 'checkbox'
                        answer{i} = val==get(h,'Max');
                     otherwise %{'float' 'integer'}
                        answer{i} = val;
                  end
               end
            case {'popupmenu' 'listbox' 'slider'}
               val = get(h,'Value');
               
               if ~isempty(val)
                  switch fmt.format
                     case 'text' % 'popupmenu' 'listbox'
                        str = get(h,'String');
                        answer{i} = str(val);
                     case 'logical' % 'checkbox'
                        answer{i} = val==get(h,'Max');
                     otherwise %{'float' 'integer'}
                        answer{i} = val;
                  end
               end
            case 'edit'
               str = get(h,'String');
               switch fmt.format
                  case {'float','integer'}
                     if ~isempty(str)
                        answer{i} = str2double(str);
                     else
                        answer{i} = [];
                     end
                  case 'date'
                     if isempty(str)
                        answer{i} = [];  % Return an empty date vector if no date was entered
                     else
                        answer{i} = datevec(str, fmt.limits);
                     end
                  case 'file'
                     if diff(fmt.limits)>1
                        answer{i} = cellstr(str);
                     else
                        answer{i} = str;
                     end
                  case {'vector'}
                     answer{i} = str2num(str); %#ok
                  otherwise %case {'text' 'dir'}
                     answer{i} = str;
               end
            case {'radiobutton' 'togglebutton'} % uibuttongroup
               hbtn = get(h,'SelectedObject');
               if isempty(hbtn) && strcmp(fmt.required,'on')
                  disp('required but not given')
               end
               if strcmp(fmt.format,'text')
                  if isempty(hbtn)
                     answer{i} = '';
                  else
                     answer{i} = get(hbtn,'String');
                  end
               else
                  if isempty(hbtn)
                     answer{i} = [];
                  else
                     answer{i} = find(hbtn==get(h,'UserData'));
                  end
               end
            case 'pushbutton'
               if strcmp(fmt.type,'color')
                  answer{i} = get(h,'BackgroundColor');
                  if strcmp(fmt.format,'integer')
                     answer{i} = uint8(round(answer{i}*255));
                  end
               else
                  answer{i} = get(h,'UserData');
               end
            case 'table'
               answer{i} = get(h,'Data');
         end
      end
   end

   function reqmet = checkReq(answer)
      idx = cellfun(@isempty,answer);
      reqmet = ~any(IsRequired(:)&idx(:));
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control Button Callback callback functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function doFigureKeyPress(obj, evd)
      [tf,I] = ismember(evd.Key,{'return','space','escape'});
      if ~tf, return; end % nothing special to do
      
      ReturnPress = false;
      if any(I==[1 2])
         
         % Ignore under a condition when GCO is an edit uicontrol
         % Known Potential Issue: GCO could be modified from outside
         hedit = findobj(gco,'flat','type','uicontrol','style','edit');
         
         ignore = ~isempty(hedit);
         if ignore
            % check for the conditions
            if I==1 % return
               % resume only if not currently in a multi-line edit uicontrol
               ignore = get(hedit,'Max')-get(hedit,'Min')>1;
               ReturnPress = true;
            end
            
            if ignore
               return;
            end
         end

         % equivalent to pressing OK button
         set(obj,'UserData','OK');
      elseif I==3 % Cancel
         % equivalent to pressing Cancel button
         set(obj,'UserData','Cancel');
      else
         return;
      end
      
      % set focus on OK button (to update the value of the current control)
      uicontrol(handles.btns(1));

      % if reached this far, valid key press to close the dialog
      uiresume(obj);

   end

   function doOKCancel(~, ~, isok)
      if isok
         set(gcbf,'UserData','OK');
         uiresume(gcbf);
      else
         set(gcbf,'UserData','Cancel');
         uiresume(gcbf); % cancel
      end
   end

   function doApply(~,~)
      DefaultAnswer = getAnswer(); % retrieve the current answers from the controls
      DefaultReqMet = checkReq(DefaultAnswer);
      Applied = true;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UICONTROL ButtonDownFcn callback functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function checkNumericRange(hObj,evt,k,cbfcn)
      
      fmt = Formats(k);
      errored = false;
      
      str = get(hObj,'String');
      if isempty(str)
         if fmt.required
            % show an error dialog
            h = errordlg('This parameter must be filled with a value.','Required Value','modal');
            errored = true;
         else
            Answer{k} = [];
         end
      else
         
         % convert to float/integer
         val = str2double(str);
         isint = strcmp(fmt.format, 'integer');
         if isint
            val = round(val);  % Round to the nearest integer
         end

         % validate the value
         lim = fmt.limits;
         if ~iscell(lim)
            lim = {'' lim(1) '' lim(2)};
            if lim{4}==inf
               lim([3 4]) = [];
            else
               lim{3} = '<=';
            end
            if lim{2}==-inf
               lim([1 2]) = [];
            else
               lim{1} = '>=';
            end
         end
         
         try
            validateattributes(val,{'numeric'},[{'scalar'} lim]);
            % Re-format the control's text according to the value
            set(hObj, 'String', num2str(val));
            Answer{k} = val; % store the numeric answer to revert to later
         catch MExcept
            h = errordlg(MExcept.message,'Invalid Value','modal');
            errored = true;
         end
      end
      
      
      if errored
         % set flag
         set(handles.fig,'UserData','Error');
         
         % revert back to the previous value
         set(hObj,'String',num2str(Answer{k}));
         uiwait(h);
      elseif ~isempty(cbfcn) % success
         % run custom callback function
         cbfcn(hObj,evt,handles.ctrls,k);
      end
      
   end

   function checkVector(hObj,evd,k,cbfcn)
      
      fmt = Formats(k);
      errored = false;
      
      str = get(hObj,'String');
      if isempty(str)
         if fmt.required
            % show an error dialog
            h = errordlg('This parameter must be filled with a value.','Required Value','modal');
            errored = true;
         else
            Answer{k} = [];
         end
      else
         % convert to float/integer
         val = str2num(str); %#ok
         
         lim = fmt.limits;
         
         try
            if iscell(lim)
               if any(strcmpi(lim,'row'))
                  type = {};
               else
                  type = {'column'};
               end
               
               validateattributes(val,{'numeric'},[type lim]);
            else
               N = numel(val);
               if N<lim(1) || N>lim(2) % incompatible dimension
                  error('This vector parameter must have %d to %d elements.',lim(1),lim(2));
               end
            end
            
            % Re-format the control's text according to the value
            Answer{k} = val; % force to be a column vector and store to revert to later
            set(hObj, 'String', num2str(Answer{k}));
         catch MExcept
            h = errordlg(MExcept.message,'Invalid Value','modal');
            errored = true;
         end
      end

      if errored
         % set flag
         set(handles.fig,'UserData','Error');
         
         % revert back to the previous value
         set(hObj,'String',num2str(Answer{k}));
         uiwait(h);
      elseif ~isempty(cbfcn) % success
         % run custom callback function
         cbfcn(hObj,evd,handles.ctrls,k);
      end
      
   end

   function checkDate(hObj,evd,k,cbfcn)

      format = Formats(k).limits;
      errored = false;
      
      str = get(hObj,'string');
      if isempty(str) % Avoid calling datenum() which prints a warning for empty strings
         if fmt.required
            % show an error dialog
            h = errordlg('This parameter must be filled with a value.','Required Value','modal');
            errored = true;
         else
            Answer{k} = '';
         end
      else
         try
            num = datenum(str, format);  % Check if the input matches the custom date format first
         catch
            try
               num = datenum(str);  % Check if the input matches any other supported date format
            catch
               h = errordlg(sprintf('Unsupported date format.'),'Invalid Value','modal');
               errored = true;
            end
         end
         Answer{k} = datestr(num,format);
         set(hObj,'String',Answer{k});
      end

      if errored
         % set flag
         set(handles.fig,'UserData','Error');
         
         % revert back to the previous value
         set(hObj,'String',num2str(Answer{k}));
         uiwait(h);
      elseif ~isempty(cbfcn)
         % run custom callback function
         cbfcn(hObj,evd,handles.ctrls,k);
      end
      
   end

   function openFilePrompt(hObj,evd,k,cbfcn)
      fmt = Formats(k);
      spec = fmt.items;
      opt = {};
      mode = diff(fmt.limits);
      filename = get(hObj,'String');
      if mode<=0 % uiputfile
         uifilefcn = @uiputfile;
         file = filename;
      else
         uifilefcn = @uigetfile;
         if mode>1 % multi-select
            file = get(hObj,'UserData');
            opt = {'MultiSelect','on'};
         else
            file = filename;
         end
      end
      
      % open the file prompt
      [f,p] = uifilefcn(spec,'',file,opt{:});
      if ~(isequal(f,0) && isequal(p,0)) % canceled, no change

         % store & display the data
         file = fullfile(p,f); % form full path(es) to the selected files
         set(hObj,'String',file);
         if mode>1
            set(hObj,'UserData',p);
         end
         
         % run custom callback function
         if ~isempty(cbfcn)
            cbfcn(hObj,evd,handles.ctrls,k);
         end
      end
      
      % bring focus back to the control
      uicontrol(hObj);
   end

   function openDirPrompt(hObj,evd,k,cbfcn)
   
      p = uigetdir(get(hObj,'String'));
      if ~isequal(p,0)
         % display the selected directory
         set(hObj,'String',p);
         
         % run custom callback function
         if ~isempty(cbfcn)
            cbfcn(hObj,evd,handles.ctrls,k);
         end
      end
      
      % bring focus back to the control
      uicontrol(hObj);
   end

   function openColorPrompt(hObj,evd,k,cbfcn)
      
      p = uisetcolor(get(hObj,'BackgroundColor'));
      if ~isempty(p)
         % display the selected color
         set(hObj,'BackgroundColor',p);
         
         % run custom callback function
         if ~isempty(cbfcn)
            cbfcn(hObj,evd,handles.ctrls,k);
         end
      end
      
      % bring focus back to the control
      uicontrol(hObj);
   end

   function forceInteger(hObj,evd,k,cbfcn)
      
      % make sure the rounded value is within the allowed range
      Answer{k} = max(get(hObj,'Min'),min(get(hObj,'Max'),round(get(hObj,'Value'))));
      
      % set the rounded value
      set(hObj,'Value',Answer{k});
      
      % run custom callback function
      if ~isempty(cbfcn)
         cbfcn(hObj,evd,handles.ctrls,k);
      end
   end
end

function S = selectivecell2struct(C,fields)
idx = ~cellfun(@isempty,fields);
if any(idx)
   idx = find(idx)';
   S = cell2struct(C(idx,:),fields(idx),1);
else
   S = C;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKPROMPT :: Check Prompt input is valid & fill default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Prompt,FieldNames,Units,TooltipStrings,err] = checkprompt(Prompt)

% default configuration
FieldNames = {}; % answer in a cell
Units = {}; % no units
TooltipStrings = {};

% standard error
err = {'inputsdlg:InvalidInput','Prompt must be a cell string with up to four columns.'};

if isempty(Prompt), Prompt = {'Input:'};
elseif ~iscell(Prompt), Prompt = cellstr(Prompt);
end

[nrow,ncol] = size(Prompt);

% prompt given in a row -> transpose
if ncol>4
   if nrow<4, Prompt = Prompt.'; [nrow,ncol] = size(Prompt);
   else return; % too many columns given
   end
end

% struct fields defined
if ncol>1
   idx = cellfun(@isempty,Prompt(:,2));
   FieldNames = Prompt(:,2);
   FieldNames(idx) = {''}; % make sure it is empty cellstr
   
   idx(:) = ~idx;
   if numel(unique(FieldNames(idx)))~=sum(idx)
      err{2} = 'Duplicate struct field name found.';
      return;
   end
else
   FieldNames = repmat({''},nrow,1);
end

% unit labels defined
if ncol>2, Units = Prompt(:,3);
else       Units = repmat({''},nrow,1);
end

% tooltip strings defined
if ncol>3, TooltipStrings = Prompt(:,4);
else       TooltipStrings = repmat({''},nrow,1);
end

% return only the labels in Prompt argument
Prompt(:,2:end) = [];

err = {}; % all cleared

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKFORMATS :: Check Formats input is valid & fill default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Formats,err] = checkformats(Formats,NumQuest,fields)

err = {};

if isempty(Formats) % if Formats not defined, use the first entry (edit:text)
   Formats = fields(1:2,:);
   Formats(cellfun(@(c)iscell(c)&&isempty(c),Formats)) = {{{}}};
   Formats = repmat(struct(Formats{:}),NumQuest,1);
end

fnames = lower(fields(1,:));
fields(1,:) = [];
nfields = numel(fnames); % sans the first row

% backward compatibility (NumLines)
if isnumeric(Formats)
   [rw,cl]=size(Formats);
   ok = rw==1;
   if ok
      OneVect = ones(NumQuest,1);
      if cl == 2, NumLines=Formats(OneVect,:);
      elseif cl == 1, NumLines=Formats(OneVect);
      elseif cl == NumQuest, NumLines = Formats';
      else ok = false;
      end
   end
   if rw == NumQuest && any(cl == [1 2]), NumLines = Formats;
   elseif ~ok
      err = {'MATLAB:inputdlg:IncorrectSize', 'NumLines size is incorrect.'};
      return;
   end
   
   % set to default edit control (column stacked)
   fields(3:end,:) = []; % all to be edit boxes
   Formats = repmat(struct(fields{:}),NumQuest,1);
   
   % set limits according to NumLines(:,1)
   numlines = mat2cell([zeros(NumQuest,1) NumLines(:,1)],ones(NumQuest,1),2);
   [Formats.limits] = deal(numlines{:});
   
   % sets the width to be 10*NumLines(:,2)
   if (size(NumLines,2) == 2)
      sizes = mat2cell([zeros(NumQuest,1) NumLines(:,2)],ones(NumQuest,1),2);
      [Formats.size] = deal(sizes{:});
   end
   
   return;
elseif ischar(Formats) || iscellstr(Formats) % given type
   if ischar(Formats), Formats = cellstr(Formats); end
   Formats = cell2struct(Formats,'type',3);
elseif ~isstruct(Formats)
   err = {'inputsdlg:InvalidInput','FORMATS must be an array of structure.'};
   return
end

% Dialog grid dimension
fdims = size(Formats);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If span field is given, fill Format struct accordingly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Formats,'span')
   
   idx = ~arrayfun(@(s)isempty(s.span),Formats); 
   if ~arrayfun(@(s,sz)isnumeric(s.span) && ...
         (numel(s.span)==2 && all(s.span==floor(s.span) & s.span>0)),Formats(idx))
      err = {'inputsdlg:InvalidInput','FORMATS.span must be 2-element vectors of positive integers.'};
      return;
   end
   if ~all(arrayfun(@(s)all(s.span<=fdims),Formats(idx)))
      err = {'inputsdlg:InvalidInput','FORMATS.span extends the control out of the grid size specified by Formats matrix.'};
      return;
   end
   
   [i,j] = find(arrayfun(@(s)any(size(s.span)>[1 1]),Formats));
   for ind = [i j].' % for each spanning control
      span = Formats(ind(1),ind(2)).span;
      % extend from left on the top-most row
      for col = ind(2)+(1:span(2)-1)
         Formats(ind(1),col).type = 'none';
         Formats(ind(1),col).limits = [0 1]; % extend from left
      end
      % extend from above for all other rows
      for row = ind(1)+(1:span(1)-1)
         for col = ind(2)+(0:span(2)-1)
            Formats(row,col).type = 'none';
            Formats(row,col).limits = [1 0]; % extend from above
         end
      end
   end
   Formats = rmfield(Formats,'span');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert Formats to cell
[~,I] = ismember(lower(fieldnames(Formats)),fnames);
if any(I==0)
   err = {'inputsdlg:InvalidFormatsField','FORMATS contains invalid field name(s).'};
   return
end
fvals = cell([prod(fdims) nfields]);
fvals(:,I) = struct2cell(Formats(:)).';
Nvals = numel(Formats);

% Mark implicitly unused cells
Iempty = cellfun(@isempty,fvals); % cell with empty format
Iunk = all(Iempty,2); % if all type/format/style are all empty, cell format is completely unknown
Inone = strcmp(fvals(:,1),'none'); % explictly specified empty cell
Nmissing = NumQuest - (Nvals-sum(Iunk)-sum(Inone)); % number of unformated prompts
if Nmissing<0
   err = {'inputsdlg:InvalidInput',sprintf('%s\n%s',...
      'FORMATS must have matching number of elements to PROMPT (exluding ''none'' type).',...
      'If .span field is used, also check for overlapping controls.')};
   return
end
if sum(Iunk)-Nmissing>0 % if extra unknown cells exist, set them to 'none'
   Imorenone = find(Iunk,sum(Iunk)-Nmissing,'last');
   fvals(Imorenone,1) = {'none'}; % set the rest of empty entries to none (spacer)
end
Iempty(:,1:2) = true; % to always copy the type & format columns

noformat = cellfun(@isempty,fields(:,2));

% avaiable styles for each Format type
styles.text = {'text'};
styles.check = {'checkbox'};
styles.edit = {'edit'};
styles.list = {'listbox','popupmenu','radiobutton','togglebutton'};
styles.range = {'slider'};
styles.table = {'table'};
styles.color = {'pushbutton'};
styles.button = {'pushbutton'};

for n = 1:Nvals % for each format entry
   [type,format,style] = deal(fvals{n,1:3});
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Check type/format/style fields (Columns 1-2) 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if there are any field without type defined, set to default
   if isempty(type) % type not specified
      Ifield = true(size(fields,1),1); % can be any field
   else % type specified
      Ifield = strcmpi(type,fields(:,1));
      if strcmpi(type,'table')
         format = []; % ignore format field
         Iempty(n,2) = false; % do not copy format field
      end
   end
   
   if isempty(format) % format specified (ignore if style=table)
      if isempty(type) && ~isempty(style) % only style specified
         Ifield(:) = Ifield & strcmpi(fields(:,3),style);
      end
   else
      Ifield(:) = Ifield & (noformat|strcmpi(fields(:,2),format));
   end

   % grab the first match
   Ifield = find(Ifield,1,'first'); 
   if isempty(Ifield)
      [i,j] = ind2sub(fdims,n);
      err = {'inputsdlg:InvalidInput',invalidformat_errormessage(i,j,type,format,fields(:,1:2))};
      return;
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Set all empty fields to type's defaults
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fvals(n,Iempty(n,:)) = fields(Ifield,Iempty(n,:));
   [type,format] = deal(fvals{n,1:2});
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check style (Column 3)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~Iempty(n,3)
      try
         fvals{n,3} = validatestring(style,styles.(type));
      catch
         err = {'inputsdlg:InvalidInput','Invalid FORMATS.style for ''range'' type must be ''slider''.'};
         return
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check items (Column 4)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % slider style: items specifies its SliderStep property
   if ~Iempty(n,4)
      items = fvals{n,4};
      if strcmp(type,'range') % type = range
         if ~any(cellfun(@(c)isnumeric(c)&&numel(c)==2&&all(c>0&c<1),items))
            err = {'inputsdlg:InvalidInput','FORMATS.items for ''range'' type must be 2 element vector with values between 0 and 1.'};
            return
         end
      elseif strcmp(type,'list')
         % check items - convert if string array or numeric array given
         if ischar(items)
            fvals{n,4} = cellstr(items);
         elseif isnumeric(items)
            fvals{n,4} = num2cell(items);
         elseif ~iscellstr(items)
            err = {'inputsdlg:InvalidInput','FORMATS.items for ''list'' type must be either a cell of strings or of numbers.'};
            return
         end
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check limits (Column 5)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % edit::date - limits specifies the date display format
   if ~Iempty(n,5)
      lims = fvals{n,5};
      if strcmp(type,'check') && strcmp(format,'text')
         if ~(iscellstr(lims) && numel(lims)==2)
            err = {'inputsdlg:InvalidInput','FORMATS.limits must be given as a two-element cellstring vector for ''check:text'' control.'};
         end
      elseif strcmp(type,'edit') && strcmp(format,'date')
         if ischar(lims) % date format string given
            try
               datestr(1,lims);
            catch
               err = {'inputsdlg:InvalidInput', 'Invalid free-form format string in FORMATS.limits for ''date'' control.'};
               return;
            end
         elseif ~(isnumeric(lims) && isscalar(lims) && any(lims==[0 1 2 6 13 14 15 16 23]))
            err = {'inputsdlg:InvalidInput','FORMATS.limits for ''edit::date'' format must be one of 0,1,2,6,13,14,15,16,23.'};
            return;
         end
      elseif strcmp(type,'table') % table - limits specifies the Column widths
         if ~iscell(lims)
            err = {'inputsdlg:InvalidInput','FORMATS.limits for ''table'' type must be given as a cell vector.'};
            return;
         end
      elseif ~((isnumeric(lims) && numel(lims)==2)...
            || (iscell(lims) && strcmp(type,'edit') && any(strcmp(format,{'float','integer','vector'}))))
         err = {'inputsdlg:InvalidInput','FORMATS.limits must be given as a two-element vector.'};
         return;
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check size (Column 6)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~Iempty(n,6)
      sz = fvals{n,6};
      if ~(isnumeric(sz) && any(numel(sz)==[1 2]) && ~any(isnan(sz)))
         err = {'inputsdlg:InvalidInput','FORMATS.size must be 1 or 2 element non-NaN vector.'};
         return
      end
      
      % if only scalar value given, set 0 as the second element
      if isscalar(sz)
         fvals{n,6}(2) = 0;
      end   
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check enable (Column 7)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~Iempty(n,7)
      try
         fvals{n,7} = validatestring(fvals{n,7},{'on','inactive','off'});
      catch
         err = {'inputsdlg:InvalidInput','FORMATS.enable must be one of {''on'',''inactive'',''off''}.'};
         return;
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check required (Column 8)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~Iempty(n,8) 
      try
         fvals{n,8} = validatestring(fvals{n,8},{'on','off'});
      catch
         err = {'inputsdlg:InvalidInput','FORMATS.required must be ''on'' or ''off''.'};
         return;
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check callback (Column 9)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~Iempty(n,9)
      cb = fvals{n,9};
      
      if isstruct(cb) % given as struct, make sure all its elements are function handles
         if all(structfun(@(cb)isa(cb,'function_handle'),cb))
            err = {'inputsdlg:InvalidInput','FORMATS.callback must be given as a function handle or a struct containing function handles.'};
            return;
         end
      elseif isa(cb,'function_handle') % convert to the struct form
         if strcmp(type,'table')
            fvals{n,9} = struct('CellEditCallback',cb);
         elseif strcmp(type,'list') && any(strcmp(style,{'radiobutton','togglebutton'})) % uibuttongroup
            fvals{n,9} = struct('SelectionChangeFcn',cb);
         elseif strcmp(style,'edit') && (strcmp(type,'color') || any(strcmp(format,{'file','dir','font'}))) % inactive edit uicontrol
            fvals{n,9} = struct('ButtonDownFcn',cb);
         else
            fvals{n,9} = struct('Callback',cb);
         end
      else
         err = {'inputsdlg:InvalidInput','FORMATS.callback must be given as a function handle or a struct containing function handles.'};
         return;
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check prompt label location (Column 10)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~Iempty(n,10)
      try
         fvals{n,10} = validatestring(fvals{n,10},{'lefttop','leftmiddle','leftbottom','topleft','topcenter','topright'});
      catch
         err = {'inputsdlg:InvalidInput','FORMATS.labelloc must be ''lefttop'', ''leftmiddle'', ''leftbottom'', ''topleft'', ''topcenter'', or ''topright''.'};
         return;
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check units label location (Column 11)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~Iempty(n,11)
      try
         fvals{n,11} = validatestring(fvals{n,11},{'righttop','rightmiddle','rightbottom','bottomleft','bottomcenter','bottomright'});
      catch
         err = {'inputsdlg:InvalidInput','FORMATS.unitsloc must be ''righttop'', ''rightmiddle'', ''rightbottom'', ''bottomleft'', ''bottomcenter'', or ''bottomright''.'};
         return;
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % check label margins (Column 12)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ~Iempty(n,12)
      v = fvals{n,12};
      if ~(isnumeric(v) && any(numel(v)==[1 2]) && ~any(isinf(v)|v<0))
         err = {'inputsdlg:InvalidInput','FORMATS.margin must be 1 or 2 element positive vector.'};
         return
      end
      if isscalar(v) % if scalar value given, use the same value for both margins
         fvals{n,12}(2) = fvals{n,12};
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gather back as Formats struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Formats = reshape(cell2struct(fvals,fnames,2),fdims);

end

function err = invalidformat_errormessage(i,j,type,format,fields)
% return error message for specifying invalid cell format

str = sprintf('Specified Format(%d,%d).',i,j);

if isempty(type) % type not specified
   Ifield = true(size(fields,1),1); % can be any field
else % type specified
   Ifield = strcmpi(type,fields(:,1));
   if ~any(Ifield) % invalid type
      types = unique(fields(:,1));
      err = sprintf('%stype field value (''%s'') is unsupported.\n',str,type);
      err = sprintf('%s\nSupported types are: %s',err,types{1});
      for n = 2:numel(types)-1
         err = sprintf('%s, %s',err,types{n});
      end
      if numel(types)>1
         err = sprintf('%s, and %s.',err,types{end});
      else
         err = sprintf('%s.',err);
      end
      return;
   end
end

I = strcmpi(fields(:,2),format);
if ~any(Ifield&I)
   formats = fields(Ifield,2);
   str = sprintf('%sformat field value (''%s'') is unsupported',str,format);
   if isempty(type)
      err = sprintf('%s.',str);
   else
      err = sprintf('%s for %s cell type.',str,type);
   end
   err = sprintf('%s\nSupported formats are: %s',err,formats{1});
   for n = 2:numel(formats)-1
      err = sprintf('%s, %s',err,formats{n});
   end
   if numel(formats)>1
      err = sprintf('%s, and %s.',err,formats{end});
   else
      err = sprintf('%s.',err);
   end
   return;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKDEFAULTS :: Check the specified default values are compatible
%%% with Formats and if one not given fill in an initial value
function [DefAns,TileMode,err] = checkdefaults(DefAns,Formats,FieldNames)
% AnsStr: 

% set the TileMode
TileMode = [1 1]; % default: no tiling
dims = size(Formats); % grid dimension
if any(dims==1) % To tile, formats must be given as a vector
   tiledim = (dims(2)==1) + 1; % if 1, tile rows; if 2, tile columns
   ansdims = size(DefAns);
   tf = ansdims>1;
   if isstruct(DefAns) && any(tf) % if vector, tile
      TileMode(tiledim) = numel(DefAns);
   elseif iscell(DefAns) && all(tf) % if matrix, tile
      TileMode(tiledim) = ansdims(2);
   end
end

if all(TileMode==1) % no tiling
   DefAns = DefAns(:);
end

% reshape to a "row" vector & trim Formats to only include relevant entries (non-'none' types)
Formats = Formats'; % go through row first
Formats = Formats(~strcmp('none',{Formats.type}));
len = numel(Formats); % if tiled, expects len*prod(TileMode) DefAns

if isempty(DefAns) % if DefAns not given
   DefAns = cell(len,1); % will set DefAns to default values
elseif isstruct(DefAns)
   if isempty(FieldNames) % FieldNames must be given via PROMPT
      err = {'inputsdlg:InvalidInput','DEFAULTANSWER is given in as a struct but its field names are not specified in PROMPT (the 2nd column).'};
      return;
   end
   
   % Convert Struct to cell according to FieldNames
   Inotempty = find(~cellfun(@isempty,FieldNames)); % ignore prompt w/o FieldName
   [~,I] = ismember(fieldnames(DefAns),FieldNames(Inotempty));
   if any(I==0)
      err = {'inputsdlg:InvalidFormatsField','DEFAULTANSWER contains invalid field name(s).'};
      return
   end
   DefStr = DefAns; % save the original input
   DefAns = cell(len,prod(TileMode)); % len==numel(FieldNames) is guaranteed
   DefAns(Inotempty(I),:) = struct2cell(DefStr(:));
   
elseif ~iscell(DefAns)
   err = {'inputsdlg:InvalidInput','Default answer must be given in a cell array or a structure.'};
   return;
elseif length(DefAns)~=len
   err = {'inputsdlg:InvalidInput','Default answer cell dimension disagrees with the number of prompt'};
   return;
end

% go through each default values
for k = 1:len
   
   fmt = Formats(k);
   
   % if any of the tiled element is given empty & its answer is required
   Iempty = cellfun(@isempty,DefAns(k,:));
   if any(Iempty) && strcmp(Formats(k).required,'on') 

      switch fmt.type
         case 'check' % off
            if strcmp(fmt.format,'logical')
               val = false;
            else
               val = fmt.limits(1);
            end
         case 'edit'
            switch fmt.format
               case {'float','integer'}
                  liminf = isinf(fmt.limits);
                  if all(liminf) % both limits inf
                     val = 0;
                  elseif any(liminf) % 1 limit inf
                     val = fmt.limits(~liminf);
                  else % neither is inf
                     val = round(mean(fmt.limits));
                  end
               case 'vector'
                  val = zeros(1,fmt.limits(1));
               otherwise %{'text','date','file','dir'}
                  val = '';
            end
         case 'list' % first item
            val = 1;
            if strcmp(fmt.format,'text')
               val = fmt.items{find(~cellfun(@isempty,fmt.items),1)};
            end
         case 'range' % middle value
            val = mean(fmt.limits);
         case 'color'
            val = [1 1 1];
         case 'table'
            val = {};
         otherwise
            val = [];
      end
      
      % set the default value to all empty controls on all tiles
      [DefAns{k,Iempty}] = deal(val);
   end

   % for all entries that are not empty check the validity of given values
   Ifilled = ~Iempty;
   if ~strcmp(fmt.type,'table') && any(Ifilled) 

      vals = DefAns(k,Ifilled); % given default values
      
      switch fmt.format
         case 'text'
            % must be char or cellsctr
            Jischar = cellfun(@ischar,vals);
            Jiscellstr = cellfun(@iscellstr,vals);
            if ~all(Jischar|Jiscellstr)
               err = {'inputsdlg:InvalidInput','Default text data format must be char.'};
               return;
            end
            
            % for the list type, value must be one of the allowed item
            if strcmp(fmt.type,'list')
               isUiControl = ismember(fmt.style,{'listbox', 'popupmenu'});
               if any(Jischar)
                  try
                     vals(:) = cellfun(@(v)validatestring(v,fmt.items),vals,'UniformOutput',false);
                  catch ME
                     err = {'inputsdlg:InvalidInput',ME.message};
                     return
                  end
                  if isUiControl
                     [~,idx] = ismember(vals,fmt.items);
                     % convert to integer format
                     vals(Jischar) = num2cell(idx);
                  end
               end
               if any(Jiscellstr)
                  Jiscellstr = find(Jiscellstr);
                  for j = Jiscellstr(:).'
                     [tf,idx] = ismember(vals{j},fmt.items);
                     if ~all(tf)
                        err = {'inputsdlg:InvalidInput','Default list item is not valid.'};
                        return;
                     end
                     % convert to integer format
                     if isUiControl
                        vals{j} = idx;
                     end
                  end
               end
                  
               DefAns(k,Ifilled) = vals;
            end
            
         case 'date' % given as a date vector
            % must be a date string or date number or date vector
            if ~all(cellfun(@(v)ischar(v)||isnumeric(v),vals))
               err = {'inputsdlg:InvalidInput','Default date date must be a valid datenum, datestr, or datevec.'};
               return;
            end
            % store the date as a string
            try
               DefAns(k,Ifilled) = cellfun(@(v)datestr(v, fmt.limits),vals,'UniformOutput',false);
            catch % not a valid date input given
               err = {'inputsdlg:InvalidInput','Default date date must be a valid datenum, datestr, or datevec.'};
               return;
            end
         case 'float'
            if strcmp(fmt.type,'color') % must be a valid RGB tuple
               attr = {'numel',3,'nonnegative','<=',1};
               varname = 'Default color data';
            else % for all other types, value must be scalar
               if iscell(fmt.limits)
                  attr = [{'scalar'} fmt.limits];
               else
                  attr = {'scalar', '>=',fmt.limits(1),'<=',fmt.limits(2)};
               end
               varname = 'Default float data';
            end
            
            try
               cellfun(@(v)validateattributes(v,{'numeric'},attr,mfilename,varname),vals);
            catch ME
               err = {'inputsdlg:InvalidInput',ME.message};
               return;
            end
            
         case 'integer' % can be multi-select if type=list
            switch fmt.type
               case 'list'
                  isUiBtnGrp = ismember(fmt.style,{'togglebutton', 'radiobutton'});
                  msel = strcmp(fmt.style,'listbox') && diff(fmt.limits)>1;
                  
                  % must be a valid index to items and
                  % if multiple-selection is not enabled, must be scalar
                  attr = {'integer','positive','<=',numel(fmt.items)};
                  if ~msel
                     attr{end+1} = 'scalar'; %#ok
                  end
                  
                  try
                     cellfun(@(v)validateattributes(v,{'numeric'},attr,mfilename,...
                        'Default list index data'),vals);
                  catch ME
                     err = {'inputsdlg:InvalidInput',ME.message};
                     return;
                  end
                  
                  if msel % if multiple-selection enabled, make sure values are unique
                     DefAns(k,Ifilled) = cellfun(@unique,vals,'UniformOutput',false);
                  elseif isUiBtnGrp
                     DefAns(k,Ifilled) = vals;%fmt.items(vals);
                  end
                  
               case 'color' % must be a RGB tuple
                  try
                     cellfun(@(v)validateattributes(v,{'numeric'},...
                        {'numel',3,'integer','nonnegative','<=',255},mfilename,...
                        'Default color data'),vals);
                  catch ME
                     err = {'inputsdlg:InvalidInput',ME.message};
                     return;
                  end
                  
                  % convert to float representation
                  DefAns(k,Ifilled) = cellfun(@(v)double(v)/255,vals,'UniformOutput',false);
               case 'check' % must be a scalar and one of limits
                  if ~all(cellfun(@(v)isscalar(v) && any(v==fmt.limits)))
                     err = {'inputsdlg:InvalidInput','Default integer check data must be a scalar value matching Format.Limits'};
                     return;
                  end
               otherwise %must be a scalar
                  % limits specifies the value range
                  if iscell(fmt.limits)
                     attr = [{'scalar','integer'} fmt.limits];
                  else
                     attr = {'scalar','integer','>=' fmt.limits(1),'<=',fmt.limits(2)};
                  end
                  try
                     cellfun(@(v)validateattributes(v,{'numeric'},attr,mfilename,'Defaultinteger data'),vals);
                  catch ME
                     err = {'inputsdlg:InvalidInput',ME.message};
                     return;
                  end
            end
            
         case 'logical'
            
            if ~all(cellfun(@(v)islogical(v)&&isscalar(v),vals))
               err = {'inputsdlg:InvalidInput','Default logical data must be of logical or numeric scalar.'};
               return;
            end
            
            % convert to integer format
            DefAns(k,Ifilled) = cellfun(@(v)fmt.limits(v+1),vals,'UniformOutput',false);
            
         case 'file'
            % must be char or cellstring
            if ~all(cellfun(@(v)ischar(v)||iscellstr(v),vals))
               err = {'inputsdlg:InvalidInput','Default file data must be given as char or cellstr.'};
               return;
            end
            
            % make everything cellstr
            vals = cellfun(@cellstr,vals,'UniformOutput',false);
            
            dlim = diff(fmt.limits);

            if dlim<=1 && ~all(cellfun(@isscalar,vals)) % single-file control
               err = {'inputsdlg:InvalidInput','Multiple default files are given for single-file control.'};
               return;
            end
            
            if dlim>=0 % for uigetfile, either directory name or existing file name
               if ~any(cellfun(@(v)all(cellfun(@(file)exist(file,'file'),v)),vals))
                  err = {'inputsdlg:InvalidInput','Default file for uigetfile must exist on the computer.'};
                  return;
               end
               
               % resolve full file name, if only name is given for a file
               % on matlab path
               changed = false;
               for n = 1:prod(TileMode)
                  files = cellfun(@which,vals{n},'UniformOutput',false);
                  I = ~cellfun(@isempty,files);
                  changed = changed || any(I);
                  vals{n}(I) = files(I);
                  
                  % multi-select files must be on a same directory
                  dirs = cellfun(@fileparts,vals{n},'UniformOutput',false);
                  if ispc
                     dirs = lower(dirs);
                  end
                  if numel(unique(dirs))>1
                     err = {'inputsdlg:InvalidInput','Default files for multi-select control must be from a same directory.'};
                     return;
                  end
               end
               
               if changed
                  DefAns(k,Ifilled) = vals;
               end
            end
         case 'dir'
            if ~all(cellfun(@(v)ischar(v)&&isrow(v)&&isdir(v),vals)) % directory must exist
               err = {'inputsdlg:InvalidInput','Default dir must be a valid path.'};
               return;
            end
         case 'vector'
            iscol = true;
            if iscell(fmt.limits)
               if any(strcmpi(fmt.limits,'row'))
                  iscol = false;
                  attr = fmt.limits;
               else
                  attr = [{'column'},fmt.limits];
                  vals{1} = vals{1}(:);
               end
               try
                  cellfun(@(v)validateattributes(v,{'numeric'},attr,mfilename,'Default vector data'),vals);
               catch ME
                  err = {'inputsdlg:InvalidInput',ME.message};
                  return;
               end
            else
               nel = numel(vals{1});
               if any(cellfun(@(v)~isnumeric(v) || nel<fmt.limits(1) || nel>fmt.limits(2),vals))
                  err = {'inputsdlg:InvalidInput','Default vector data must be numeric and the number of elements must be within the limit.'};
                  return;
               end
            end
            if iscol
               idx = cellfun(@(v)~iscolumn(v),vals);
               DefAns(k,Ifilled) = cellfun(@(v)v(:),vals(idx),'UniformOutput',false);
            end
      end
   end
end

DefAns = DefAns(:);
err = {}; % all good
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Options,FormatFields,err] = checkoptions(UserOptions)

err = {'inputsdlg:InvalidInput',''};

% default options
Fields = {
   'Resize'             'off'
   'WindowStyle'        'normal'
   'Interpreter'        'tex'
   'DialogSizeMinimum'	[inf inf] % automatically set
   'CancelGiven'        false
   'CancelButton'       'on'
   'ApplyGiven'         false
   'ApplyButton'        'off'
   'ButtonNames'        {{'OK','Cancel','Apply'}}
   'Sep'                10
   'AlignControls'      'off'
   'FontSize'           get(0,'DefaultUicontrolFontSize')
   'CreateFcn'          {{}}
   'DeleteFcn'          {{}}
   }.';

% default formats, given type & format
FormatFields = [
   {'type'   'format'  'style'      'items' 'limits'   'size' 'enable' 'required' 'callback' 'labelloc' 'unitsloc' 'margin'}
   {'edit'   'text'    'edit'       {}      [0 1]      [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]} % default if Formats or Formats.type not given
   {'edit'   'integer' 'edit'       {}      [-inf inf] [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'edit'   'float'   'edit'       {}      [-inf inf] [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'edit'   'vector'  'edit'       {}      [0 inf]    [0 0]  'on'     'off'      struct([]) 'lefttop'  'righttop' [3 3]}
   {'edit'   'date'    'edit'       {}      2          [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'edit'   'file'  'edit'  {'*.*' 'All Files'} [0 1] [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'edit'   'dir'     'edit'       {}      [0 1]      [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]} % default if Formats or Formats.type not given
   {'check'  'logical' 'checkbox'   {}      [0 1]      [0 0]  'on'     'on'       struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'check'  'integer' 'checkbox'   {}      [0 1]      [0 0]  'on'     'on'       struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'check'  'text'    'checkbox'   {}      {'off','on'} [0 0] 'on'    'on'       struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'list'   'integer' 'popupmenu'  {}      [0 1]      [0 0]  'on'     'on'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'list'   'text'    'popupmenu'  {}      [0 1]      [0 0]  'on'     'on'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'range'  'float'   'slider'     {}      [0 1]      [0 0]  'on'     'on'       struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'range'  'integer' 'slider'     {}      [0 255]    [0 0]  'on'     'on'       struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'color'  'float'   'pushbutton' {}      [0 1]      [65 0]  'on'     'on'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'color'  'integer' 'pushbutton' {}      [0 255]    [65 0]  'on'     'on'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'button' ''        'pushbutton' {}      [0 1]      [0 0]  'on'     'off'      struct([]) 'leftmiddle'  'righttop' [3 3]}
   {'text'   ''        'text'       {}      [0 1]      [0 0]  'on'     'off'      struct([]) 'lefttop'  'righttop' [3 3]}
   {'none'   ''        ''           {}      [0 0]      [0 0]  'on'     'off'      struct([]) 'lefttop'  'righttop' [3 3]} % default if Formats.type empty
   {'table'  {}        'table'      {}      {'auto'}   [0 0]  'on'     'off'      struct([]) 'lefttop'  'righttop' [3 3]}];

Options = struct(Fields{:});

if isempty(UserOptions) % no option specified, use default
   err = {};
   return;
elseif ischar(UserOptions) && any(strcmpi(UserOptions,{'on','off'}))
   warning off MATLAB:warn_r14_stucture_assignment
   UserOptions.Resize = UserOptions;
   warning on MATLAB:warn_r14_stucture_assignment
elseif ~isstruct(UserOptions)
   err{2} = 'Options must be ''on'', ''off'', or a struct.';
   return;
end

if numel(UserOptions)~=1
   err{2} = 'Options struct must be a scalar.';
   return;
end

% check if User Resize Option is given as on/off string

% to-do: Separate overall options to format options
% optfnames = fieldnames(UserOptions);
% optnames = regexpi(optfnames,'default(\S+)','tokens');
% Iopt = ~cellfun(@isempty,optnames);
% optnames = lower(cellfun(@(c)c{1},[optnames{Iopt}],'UniformOutput',false));
% fmtnames = lower(FormatFields(1,:));
% nfmtfields = numel(fmtnames); % sans the first row

% remove all the format options
% Options = rmfield(Options,optfnames(Iopt)));

% start with 'type'-'format' combo for the default control type
% [tf,I] = ismember(optnames,fnames);


% check UserOptions struct & update Options fields
for fname_cstr = fieldnames(UserOptions)' % for each user option field
   
   fname = char(fname_cstr); % use plain char string (not cellstr)
   val = UserOptions.(fname);
   
   % make sure string value is given as a cellstr
   if ischar(val), val = cellstr(val); end
   
   % if field not filled, use default value
   if isempty(UserOptions.(fname)), continue; end
      
   switch lower(fname)
      case 'resize'
         if numel(val)~=1 || ~any(strcmpi(val,{'on','off'}))
            err{2} = 'Resize option must be ''on'' or ''off''.';
            return;
         end
         Options.Resize = char(val);
      case 'windowstyle'
         if numel(val)~=1 || ~any(strcmpi(val,{'normal','modal','docked'}))
            err{2} = 'WindowStyle option must be ''normal'' or ''modal''.';
            return;
         end
         Options.WindowStyle = char(val);
      case 'interpreter'
         if numel(val)~=1 || ~any(strcmpi(val,{'latex','tex','none'}))
            err{2} = 'Interpreter option must be ''latex'', ''tex'', or ''none''.';
            return;
         end
         Options.Interpreter = char(val);
      case 'cancelbutton'
         if numel(val)~=1 || ~any(strcmpi(val,{'on','off'}))
            err{2} = 'CancelButton option must be ''on'' or ''off''.';
            return;
         end
         Options.CancelGiven = true;
         Options.CancelButton = char(val);
      case 'applybutton'
         if numel(val)~=1 || ~any(strcmpi(val,{'on','off'}))
            err{2} = 'ApplyButton option must be ''on'' or ''off''.';
            return;
         end
         Options.ApplyGiven = true;
         Options.ApplyButton = char(val);
      case 'buttonnames'
         if ~iscellstr(val)
            err{2} = 'ButtonNames option must be of cellstr or char type.';
            return;
         end
         
         % if not all 3 button names are given, use default for unspecified
         N = numel(val);
         if (N>3)
            err{2} = 'ButtonNames option takes up to 3 button names.';
            return;
         end
         Options.ButtonNames(1:N) = val;
      case 'sep'
         if numel(val)~=1 || ~isnumeric(val) || val<0
            err{2} = 'Sep option must be non-negative scalar value.';
            return;
         end
         Options.Sep = val;
      case 'aligncontrols'
         if numel(val)~=1 || ~any(strcmpi(val,{'on','off'}))
            err{2} = 'AlignControls option must be ''on'' or ''off''.';
            return;
         end
         Options.AlignControls = char(val);
      case 'unitsmargin'
         error('UnitMargin option has been deplicated. Use Formats.margin to set individual control''s label margins.');
      case 'fontsize'
         if ~(isscalar(val) && isnumeric(val) && val>0)
            err{2} = 'FontSize option must be non-negative scalar value.';
            return;
         end
         Options.FontSize = val;
      case 'dialogsizeminimum'
         if ~(isnumeric(val) && numel(val)==2 && ~any(isnan(val)))
            err{2} = 'DialogSizeMinimum option must be 2-element vector.';
         end
         % no minimum if not positive
         idx = val<=0;
         Options.DialogSizeMinimum(idx) = inf;
      case 'createfcn'
         if ~(isempty(val) || isa(val,'function_handle'))
            err{2} = 'CreateFcn option must be a function handle.';
         end
         Options.CreateFcn = val;
      case 'deletefcn'
         if ~(isempty(val) || isa(val,'function_handle'))
            err{2} = 'DeleteFcn option must be a function handle.';
         end
         Options.DeleteFcn = val;
      otherwise
         warning('inputsdlg:InvalidOption','%s is not a valid option name.',fname);
   end
end

err = {}; % all cleared
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BUILDGUI :: Builds the dialog box and returns handles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Formats,handles,sinfo] = buildgui(Title,Prompt,Units,TooltipStrings,FieldNames,Formats,TileMode,Options)
% 1. Create handle graphic objects for all controls (embedded in uipanels)
% 2. Generate sinfo to assist object positioning in doResize

sep = Options.Sep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tile the controls
coltile = TileMode(2)>1;
if any(TileMode>1)
   
   Formats = repmat(Formats,TileMode);

   num = numel(Prompt)*prod(TileMode);
   Prompt = reshape(repmat(Prompt,TileMode),num,1);
   FieldNames = reshape(repmat(FieldNames,TileMode),num,1);
   Units = reshape(repmat(Units,TileMode),num,1);
   TooltipStrings = reshape(repmat(TooltipStrings,TileMode),num,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% determine how to utilize 'none' space
% Place all the elements at (0,0)
dim = size(Formats); % display grid dimension
free = strcmp('none',{Formats.type}); % location of empty block maybe to be occupied by neighbor entry
num = sum(~free); % number of controls
map = zeros(dim); % determine which control occupies which block(s)
order = zeros(1,num); % uicontrol placement order (for tab order)
wcell = zeros(dim); % -> widths of cells
hcell = zeros(dim); % -> heights of cells
autoheight = false(dim(1),1);
autowidth = false(1,dim(2));

n = 1;
for f = 1:prod(dim)
   if coltile % traverse column-first
      [i,j] = ind2sub(dim,f);
   else % traverse row-first
      [j,i] = ind2sub(dim([2 1]),f);
   end
   m = sub2ind(dim,i,j);
   
   if free(m)
      mode = diff(Formats(m).limits);
      [i,j] = ind2sub(dim,m);

      if mode>0 && j>1, map(m) = map(sub2ind(dim,i,j-1)); % copy from left
      elseif mode<0 && i>1, map(m) = map(sub2ind(dim,i-1,j)); % copy from above
      end % other wise, 0 (nothing occupying)
   else
      map(m) = n;
      order(n) = m;
      n = n + 1;
   end
end

% Check none-type cells, possibly acting as a placeholder
[I,J] = find(reshape(free,dim));
for n = 1:numel(I)
   sz = Formats(I(n),J(n)).size;
   if sz(1)>0 % placeholding for fixed row height
      hcell(I(n),J(n)) = sz(1);
   elseif sz(1)<0 % placeholding with minimum row height
      hcell(I(n),J(n)) = -sz(1);
      autoheight(I(n)) = true;
   end
   if sz(2)>0 % placeholding for fixed column width
      wcell(I(n),J(n)) = sz(2);
   elseif sz(2)<0 % placeholding with minimum row height
      wcell(I(n),J(n)) = -sz(2);
      autowidth(I(n)) = true;
   end
end

% remove none's from Formats and order the rest in Prompt order
Formats = Formats(order).'; % removes all none-types

FigColor=get(0,'DefaultUicontrolBackgroundcolor');

fig = dialog(           ...
   'Visible'     ,'off'   , ...
   'Name'       ,Title   , ...
   'Pointer'     ,'arrow'  , ...
   'Units'      ,'pixels'  , ...
   'UserData'     ,'Cancel'  , ...
   'Tag'       ,'Inputsdlg'   , ...
   'HandleVisibility' ,'callback' , ...
   'Color'      ,FigColor  , ...
   'WindowStyle'   ,Options.WindowStyle, ...
   'DoubleBuffer'   ,'on'    , ...
   'Resize'      ,Options.Resize    ...
   );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create controls
%%%%%%%%%%%%%%%%%%%%%
CommonInfo = {'Units'  'pixels'
   'FontSize'      Options.FontSize
   'FontWeight'   get(0,'DefaultUicontrolFontWeight');
   'HandleVisibility'  'callback'}';

props.edit = [CommonInfo {...
   'Style'      'edit';
   'HorizontalAlignment' 'left';
   'BackgroundColor' 'white'}'];

props.checkbox = [CommonInfo {...
   'Style'      'checkbox';
   'HorizontalAlignment' 'left';
   'BackgroundColor' FigColor}'];

props.popupmenu = [CommonInfo {...
   'Style'      'popupmenu';
   'HorizontalAlignment' 'left';
   'BackgroundColor' 'white'}'];

props.listbox = [CommonInfo {...
   'Style'      'listbox';
   'HorizontalAlignment' 'left';
   'BackgroundColor' 'white'}'];

props.slider = [CommonInfo {...
   'Style'      'slider';
   }'];

props.uibuttongroup = [CommonInfo {...
   'BackgroundColor' FigColor
   }'];

props.radiobutton = props.checkbox;
props.radiobutton{2,strcmp(props.checkbox(1,:),'Style')} = 'radiobutton';

props.pushbutton = [CommonInfo {...
   'Style'        'pushbutton';
   'HorizontalAlignment' 'center'}'];

props.togglebutton = props.pushbutton;
props.togglebutton{2,strcmp(props.pushbutton(1,:),'Style')} = 'togglebutton';

props.table = CommonInfo;

% Add VerticalAlignment here as it is not applicable to the above.
props.text = [CommonInfo {...
   'Style'               'text'
   'Units'               'normalized'
   'Position'            [0 0 1 1]
   'Visible'             'off'}'];

props.label = [CommonInfo {...
   'BackgroundColor'     FigColor;
   'HorizontalAlignment' 'left';
   'VerticalAlignment'   'bottom';
   'Color'        get(0,'FactoryUIControlForegroundColor');
   'Interpreter'     Options.Interpreter}'];

% For each control, place a uicontrol and an axes (axes is to enable LaTeX)
figsz = get(fig,'Position');
props.uipanel = [CommonInfo {...
   'BackgroundColor' FigColor
   'BorderType','none'
   }.'];
props.axes = [CommonInfo {...
   'Units' 'normalized'
   'Position' [0 0 1 1]
   'XLim' [0 figsz(3)]
   'YLim' [0 figsz(4)]
   'Visible' 'off'}.'];

isbtngrp = false(num,1);
istext = reshape(arrayfun(@(s)strcmp(s.type,'text'),Formats),num,1);
autosize = zeros(num,2); % 0-fixed, 1-autosize, 2-resize with window
rowspan = zeros(num,1); % # of rows control occupies
colspan = zeros(num,1); % # of columns control occupies
hPanels = zeros(num,2); % [uipanel|axes]
hCtrls = zeros(num,3); % [Prompt|Edit|Unit]
for m = 1:num % for each control
   
   % get current control's Format spec
   fmt = Formats(m);
   
   % determine the span of the control
   [i,j] = find(map==m);
   rowspan(m) = numel(unique(i));
   colspan(m) = numel(unique(j));
   
   % set autosize (if fixed, set only once)
   %    autosize(m,:) = (fmt.size<0) + (fmt.size<=0);
   autosize(m,:) = (fmt.size<=0);
   
   % if autosize, set to default size, if fixed
   if autosize(m,1) % width
      if fmt.size(1)==0
         fmt.size(1) = 1; % temp size
      else
         fmt.size(1) = -fmt.size(1);
      end
   end
   if autosize(m,2) % height
      if fmt.size(2)==0
         fmt.size(2) = 1; % temp size
      else
         fmt.size(2) = -fmt.size(2);
      end
   end
   
   % for uicontrols
   hPanels(m,1) = uipanel('Parent',fig,'Position',[0 0 fmt.size],props.uipanel{:});
   
   % for text controls & labels
   hPanels(m,2) = axes(props.axes{:},'Parent',hPanels(m,1));
   
   % always add labels even if control is not using them
   hCtrls(m,2) = text('Parent',hPanels(m,2),props.label{:},'String',Prompt{m});
   hCtrls(m,3) = text('Parent',hPanels(m,2),props.label{:},'String',Units{m});
   
   idx = strcmp(fmt.style,{'radiobutton','togglebutton'});
   isbtngrp(m) = any(idx);
   if isbtngrp(m)
      % create the UI Button Group object
      h = uibuttongroup('Parent',hPanels(m,1),props.uibuttongroup{:},...
         'Position',[0 0 fmt.size],'Tag',FieldNames{m});%,'Title',Prompt{m}
      hCtrls(m,1) = h;
      
      % Create button objects and record their extents
      dim_btns = size(fmt.items);
      kvalid = find(~cellfun(@isempty,fmt.items));
      Nvalid = numel(kvalid);
      hButtons = zeros(Nvalid,1);
      btn_w = zeros(dim_btns);
      btn_h = zeros(dim_btns);
      for n = 1:numel(kvalid)
         [i,j] = ind2sub(dim_btns,kvalid(n));
         hButtons(n) = uicontrol('Parent',h,'Style',fmt.style,props.(fmt.style){:},...
            'String',fmt.items{i,j},'Min',0,'Max',1,'UserData',kvalid(n),...
            'TooltipString',TooltipStrings{m});
         pos = get(hButtons(n),'Extent');
         btn_w(i,j) = pos(3);
         btn_h(i,j) = pos(4);
      end
      
      % set buttons sizes and extra margins to account for the button width
      if idx(1) % radiobutton
         % each column to have the same width
         margin = [20 0];
         btn_w = max(btn_w,[],1) + margin(1); % button widths with extra pixels for non-text part of control
      else % togglebutton
         % all buttons to have the same width
         margin = [12 2];
         btn_w = repmat(max(btn_w(:)) + margin(1),1,dim_btns(2)); % button widths with extra pixels for non-text part of control
      end
      btn_h = max(btn_h,[],2) + margin(2); % button heights with extra pixels for non-text part of conrol
      
      % Button positions
      btn_sep = sep*3/4;
      x0 = cumsum([0 btn_w]+btn_sep)-btn_sep*3/8;
      y0 = flipud(cumsum([0;btn_h]+btn_sep)-btn_sep*3/8);
      
      % set positions of buttons
      kvalid = find(hButtons~=0);
      for n = 1:Nvalid
         [i,j] = ind2sub(dim_btns,kvalid(n)); % i-col, j-row
         pos = [x0(j) y0(i+1) btn_w(j) btn_h(i)];
         set(hButtons(n),'Position',pos);
      end
      
      % set the size
      set(h,'Position',[0 0 x0(end) y0(1)],'UserData',hButtons);
      autosize(m,:) = 0; % no autosize
      
   elseif strcmp(fmt.style,'table') % uitable
      
      hCtrls(m,1) = uitable('Parent',hPanels(m,1), props.(fmt.style){:}, ...
         'Position',[0 0 fmt.size],'ColumnName',fmt.items,...
         'ColumnFormat',fmt.format,'ColumnWidth',fmt.limits,...
         'ColumnEditable',true,...
         'Tag',FieldNames{m});
      
   else % uicontrols
      
      % create a uipanel, embed a uicontrol, and 2 text labels
      hc = uicontrol('Parent',hPanels(m,1), 'Style',fmt.style, props.(fmt.style){:},...
         'Position',[0 0 fmt.size],'Tag',FieldNames{m},'TooltipString',TooltipStrings{m});
      hCtrls(m,1) = hc;
      
      % set min and max if not a numeric edit box
      if ~(any(strcmp(fmt.style,'edit')) && any(strcmp(fmt.format,{'float','integer','date','dir','vector'}))) ...
            && ~(strcmp(fmt.type,'check')&&strcmp(fmt.format,'text'))
         lim = fmt.limits;
         if any(isinf(lim)), lim = [0 2]; end % edit:vector
         set(hc,'Min',lim(1),'Max',lim(2));
      end
      
      % style-dependent configuration
      switch fmt.style
         case 'text' % static text (only a label)
            
            % display the label on the upper left hand corner of the display grid cell
            set(hCtrls(m,2),'Units','normalized','Position',[0 1],'VerticalAlignment','top');
            
            % create invisible dummy control
            set(hc,'FontName',get(hCtrls(m,2),'FontName'),'String',Prompt{m});
            
            % if not autowidth, go ahead and map-out
            if ~autosize(m,1)
               set(hc,'Position',[0 0 fmt.size]);
               msg = textwrap(hc,Prompt(m));
               str = sprintf('%s\n',msg{:});
               set(hc,'String',str(1:end-1));
               autosize(m,2) = 0; % only auto-height if auto-width
            end
            
            % Use Axes Text object (in order to render LaTeX text)
            set(hCtrls(m,3),'String','','Visible','off');
            
         case 'edit' % edit type
            
            % check if multi-line control
            if iscell(fmt.limits)
               dlim = 0;
            else
               dlim = round(diff(fmt.limits));
            end
            
            % set multi-line edit box
            multiline = (any(strcmp(fmt.format,{'text','file'})) && dlim>1);
            if multiline
               nrows = dlim-1;
            else
               multiline = strcmp(fmt.format,'vector');
               if multiline
                  if iscell(fmt.limits)
                     multiline = ~any(strcmpi(fmt.limits,'row')); % if row-vector, a single-line editbox suffices
                     if multiline
                        nrows = find(strcmpi(fmt.limits,'numel'),1);
                        if isempty(nrows)
                           nrows = 5;
                        else
                           nrows = min(5,fmt.limits{nrows+1});
                        end
                     end
                  else
                     nrows = max(fmt.limits(1),min(fmt.limits(2),5));
                  end
               end
            end
            if multiline
               set(hc,'Min',0,'Max',2);
            end
            
            % change alignment for numeric formats (default: center)
            if any(strcmp(fmt.format,{'float','integer'}))
               set(hc,'HorizontalAlignment','center');
            elseif strcmp(fmt.format,'vector')
               set(hc,'HorizontalAlignment','left');
            end
            
            % set auto-height (no resize allowed if single-line)
            if autosize(m,2)>0 % auto-height adjustment
               % set
               if multiline % set to have dlim lines
                  set(hc,'String',repmat(sprintf(' \n'),1,nrows-1));
               else % single-line
                  set(hc,'String',' ');
                  autosize(m,2) = 0; % no need to adjust dynamically
               end
               ext = get(hc,'Extent');
               set(hc,'String','');
               fmt.size(2) = ext(4); % set to the font height
               set(hc,'Position',[0 0 fmt.size]);
            end
            
         case 'checkbox' % no labels
            
            % Show the prompt label with the control
            set(hc,'String',Prompt{m});
            set(hCtrls(m,2),'String','','Visible','off');
            
            % Set the control size (fixed, no resize)
            pos = get(hc,'Extent'); % width reflects only label width
            pos(3) = pos(3) + 20; % pad extra for the checkbox itself
            set(hc,'Position',pos);
            autosize(m,:) = 0; % no resizing
            
         case 'popupmenu'
            
            % get the width of the widest entry
            if autosize(m,1) % auto-width
               w = 0;
               for n = 1:numel(fmt.items)
                  set(hc,'String',fmt.items{n});
                  ext = get(hc,'Extent');
                  w = max(w,ext(3));
               end
               fmt.size(1) = w + 20; % additional width for the pulldown
            end
            
            if autosize(m,2) % auto-height
               % only list 1 item and get the extent
               set(hc,'String',fmt.items{1});
               ext = get(hc,'Extent');
               fmt.size(2) = ext(4);
            end
            
            % re-set position
            if any(autosize(m,:))
               set(hc,'Position',[0 0 fmt.size]);
               autosize(m,:) = 0; % no resizing
            end
            
            % Set menu & choose the first entry
            set(hc, 'String',fmt.items);
            
         case 'listbox'
             
            % Set menu & choose the first entry
            set(hc,'String',fmt.items);
           
            if any(autosize(m,:))
               % determine the optimal size
               ext = get(hc,'Extent');
               if autosize(m,1) % auto-width
                  fmt.size(1) = ext(3) + 20;
                  if autosize(m,2)~=1 && fmt.size(2)<ext(4) % with vertical scroller
                     fmt.size(1) = fmt.size(1); % add vertical scroll bar width
                  end
                  autosize(m,1) = 0; % no resizing
               end
               if autosize(m,2) % auto-height -> set to the tallest
                  % Restrict the height if a maximum number of lines was specified
                  if fmt.limits(1) > 0 && fmt.limits(1) < numel(fmt.items)
                     set(hc,'String',fmt.items(1:fmt.limits(1)));
                     ext = get(hc,'Extent');
                     set(hc,'String',fmt.items);
                  end
                  fmt.size(2) = ext(4);
                  if fmt.limits(1)>0
                     autosize(m,2) = 0;
                  end
               end
               
               % re-set position
               set(hc,'Position',[0 0 fmt.size]);
            end
            
         case 'slider'
            
            if any(autosize(m,:))
               if autosize(m,1) && ~autosize(m,2) % vertical slider, auto-width -> fixed width
                  fmt.size(1) = 16;
                  autosize(m,1) = 0;
               elseif autosize(m,2) % auto-height -> fixed height
                  fmt.size(2) = 16;
                  autosize(m,2) = 0;
               end
               set(hc,'Position',[0 0 fmt.size]);
            end
            
            % set slider step if items field is filled
            if ~isempty(fmt.items)
               set(hc,'SliderStep',fmt.items);
            end
            
         case 'pushbutton' % button & color types

            if strcmp(fmt.type,'color') % color type
               set(hc,'String',' ') % space to set auto-height
            else % button type
               % Show the prompt label with the control
               set(hc,'String',Prompt{m});
               set(hCtrls(m,2),'String','','Visible','off');
            end
            
            if autosize(m,2) % auto-height -> fix
               ext = get(hc,'Extent');
               fmt.size(2) = ext(4) + 6;
               set(hc,'Position',[0 0 fmt.size]);
               autosize(m,2) = 0;
            end
      end
   end
end % for m = 1:num

% layout each control and labels within the panel (and axes)
labelbase = zeros(num,2,2); % lower-left corner coordinates
labelpos = zeros(num,2,2);
labelsiz = zeros(num,2,2);
labelalign = ones(num,4);
ctrlpos = zeros(num,4);
set(hCtrls(istext,2),'Units','pixels');
alignprops = cell(2,2); % {'HorizontalAlignemt','VerticalAlignment'}
for m = 1:num
   
   fmt = Formats(m);
   ext = cell2mat(get(hCtrls(m,[2 3]),'Extent'));
   if istext(m) % text type does not use uicontrol
      if autosize(m,1) % if autosize, minimum size to be the square
         ext(1,3) = ext(1,4);
      end
      pos = [0 0 0 0];
   else
      pos = get(hCtrls(m,1),'Position');
      if isbtngrp(m)
         pos(1) = pos(1) + 1;
      end
   end
   wd = [ext(1,3);pos(3);ext(2,3)];
   ht = [ext(1,4);pos(4);ext(2,4)];
   x0 = zeros(2,1); x = x0;
   y0 = zeros(2,1); y = y0;
   
   % place prompt label w.r.t. the control @ (0,0)
   tok = regexp(fmt.labelloc,'^(left|top)(.+)$','tokens','once');
   labelalign(m,1) = strcmp(tok{1},'left');
   if labelalign(m,1) % left
      % if control height is autoadjusted, make sure it is higher than label height
      if autosize(m,2)
         ht(2) = max(ht(2),ht(1));
      end
      
      alignprops{1,1} = 'left'; alignprops(1,2) = tok(2);
      labelalign(m,2) = find(strcmp(tok{2},{'bottom','middle','top'}))-1;
      x0(1) = -wd(1) - fmt.margin(1); % move label to left
      x(1) = x0(1);
      y0(1) = (ht(2)-ht(1))*labelalign(m,2)/2;
      y(1) = ht(2)*labelalign(m,2)/2;
      
   else % above
      % if control width is autoadjusted, make sure it is wider than label width
      if autosize(m,1)
         wd(2) = max(wd(2),wd(1));
      end
      
      alignprops{1,2} = 'bottom'; alignprops(1,1) = tok(2);
      labelalign(m,2) = find(strcmp(tok{2},{'left','center','right'}))-1;
      x0(1) = (wd(2)-wd(1))*labelalign(m,2)/2;
      x(1) = wd(2)*labelalign(m,2)/2;
      y0(1) = ht(2) + fmt.margin(1); % move label up
      y(1) = y0(1);
   end
   
   % units label w.r.t. the control @ (0,0)
   tok = regexp(fmt.unitsloc,'^(right|bottom)(.+)$','tokens','once');
   labelalign(m,3) = strcmp(tok{1},'right');
   if labelalign(m,3) % right
      % if control height is autoadjusted, make sure it is higher than unit height
      if autosize(m,2)
         ht(2) = max(ht(2),ht(3));
      end
      
      alignprops{2,1} = 'left'; alignprops(2,2) = tok(2);
      labelalign(m,4) = find(strcmp(tok{2},{'bottom','middle','top'}))-1;
      x0(2) = wd(2) + fmt.margin(2); % move label right
      x(2) = x0(2);
      y0(2) = (ht(2)-ht(3))*labelalign(m,4)/2;
      y(2) = ht(2)*labelalign(m,4)/2;
   else % below
      % if control width is autoadjusted, make sure it is wider than unit width
      if autosize(m,1)
         wd(2) = max(wd(2),wd(3));
      end
      
      alignprops{2,2} = 'top'; alignprops(2,1) = tok(2);
      labelalign(m,4) = find(strcmp(tok{2},{'left','center','right'}))-1;
      x0(2) = (wd(2)-wd(3))*labelalign(m,4)/2;
      x(2) = wd(2)*labelalign(m,4)/2;
      y0(2) = -ht(3)-fmt.margin(2); % move label down
      y(2) = -fmt.margin(2);
   end
   
   % translate so that all coordinates are on the first quadrant
   xmin = min(min(x0),0);
   ymin = min(min(y0),0);
   
   labelbase(m,1,:) = x0 - xmin;
   labelbase(m,2,:) = y0 - ymin;
   labelpos(m,1,:) = x - xmin;
   labelpos(m,2,:) = y - ymin;
   labelsiz(m,:,:) = ext(:,[3 4]).';
   
   ctrlpos(m,:) = [-xmin -ymin wd(2) ht(2)];
   if isbtngrp(m)
      ctrlpos(m,1) = ctrlpos(m,1) + 1;
      ctrlpos(m,3) = ctrlpos(m,3) - 1;
   end
   
   if ~istext(m)
      set(hCtrls(m,[2 3]),{'HorizontalAlignment','VerticalAlignment'},alignprops);
   end
end
set(hCtrls(istext,2),'Units','normalized');

% minimum panel size
pnpos = [zeros(num,2) max(max(labelbase+labelsiz,[],3),ctrlpos(:,[1 2])+ctrlpos(:,[3 4]))];

% if strcmpi(Options.AutoLayout,'on')
%    % minimum panel size
%    pnpos = [zeros(num,2) max(max(labelbase+labelsiz,[],3),ctrlpos(:,[1 2])+ctrlpos(:,[3 4]))];
%    map = autolayout(pnpos,Options.DesiredFigureWidth,sep);
%    dim = size(map);
% end

% Optionally align adjacent controls in each column
if strcmpi(Options.AlignControls, 'on')
   for n = 1:dim(2) % for each column
      if n==1
         notspanned = true(dim(1),1);
      else
         notspanned = diff(map(:,[n-1 n]),[],2)~=0;
      end
      idx = setdiff(unique(map(notspanned,n)),0); % must be a control & remove duplicates
      x0 = max(ctrlpos(idx,1));
      dx = x0-ctrlpos(idx,1);
      ctrlpos(idx,1) = ctrlpos(idx,1) + dx;
      labelpos(idx,1,:) = bsxfun(@plus,labelpos(idx,1,:),dx);
      labelbase(idx,1,:) = bsxfun(@plus,labelbase(idx,1,:),dx);
   end
   
   % recompute the minimum panel sizes
   pnpos = [zeros(num,2) max(max(labelbase+labelsiz,[],3),ctrlpos(:,[1 2])+ctrlpos(:,[3 4]))];
end

% position the controls w/in their respective panel, and set panel size
nottext = ~istext;
Nnottext = sum(nottext);
set(hCtrls(nottext,1),{'Position'},mat2cell(ctrlpos(nottext,:),ones(Nnottext,1),4));
set(hCtrls(nottext,2),{'Position'},mat2cell(labelpos(nottext,:,1),ones(Nnottext,1),2));
set(hCtrls(:,3),{'Position'},mat2cell(labelpos(:,:,2),ones(num,1),2));
set(hPanels(:,1),{'Position'},mat2cell(pnpos,ones(num,1),4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the maximum and minimum column widths & row heights

idx = map>0;
wctrl = (pnpos(:,3)-sep*(colspan-1))./colspan;
wcell(idx) = wctrl(map(idx));
hctrl = (pnpos(:,4)-sep*(rowspan-1))./rowspan; % heights of controls
hcell(idx) = hctrl(map(idx)); % -> heights of cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine which columns & rows to autosize. Only column or row is
% autosized for spanned controls, favoring the left-most column or
% bottom-most row

% get preferred column to autosize
% Nauto = autosize(:,1); 
Nauto = false(dim);
Nauto(idx) = autosize(map(idx),1);
Nauto = sum(Nauto,1); % number of autowidth controls in each column
[~,col] = sort(fliplr(Nauto),'descend');
col(:) = dim(2)-col+1; % preferred order of columns to be autosized

% assign which column to autosize for each control
Ictrl = find(autosize(:,1)>0); % controls with auto-sized width
for m = Ictrl.' % for each control
   [~,j] = find(map==m);
   [tf,I] = ismember(col,j);
   autowidth(j(I(find(tf,1)))) = true;
end

% get preferred row to autosize
Nauto = false(dim);
Nauto(idx) = autosize(map(idx),2);
Nauto = sum(Nauto,2); % number of autowidth controls in each column
[~,row] = sort(flipud(Nauto),'descend');
row(:) = dim(1)-row+1; % preferred order of columns to be autosized

% assign which row to autosize for each control
Ictrl = find(autosize(:,2)>0&~istext);
for m = Ictrl.'
   [i,~] = find(map==m);
   [tf,I] = ismember(row,i);
   autoheight(i(I(find(tf,1)))) = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the button panel

% All buttons on a uipanel
hBtnPanel = uipanel('Parent',fig,props.uipanel{:});

% OK Button
hBtn(1) = uicontrol(hBtnPanel, props.pushbutton{:}, ...
   'String', Options.ButtonNames{1}, 'UserData', 'OK');

editable = ~all(istext);

% Cancel Button
if (Options.CancelGiven && strcmpi(Options.CancelButton,'on')) || (~Options.CancelGiven && editable)
   hBtn(2) = uicontrol(hBtnPanel, props.pushbutton{:}, ...
      'String', Options.ButtonNames{2}, 'UserData', 'Cancel');
end

% Apply Button
if strcmpi(Options.ApplyButton,'on')
   hBtn(end+1) = uicontrol(hBtnPanel, props.pushbutton{:}, ...
      'String', Options.ButtonNames{3}, 'UserData', 'Apply');
end

% set size
offset = 25;
minwidth = 37;
Nbtns = numel(hBtn);
ext = cell2mat(get(hBtn,{'Extent'}));
btnw = max(minwidth,max(ext(:,3)))+offset;
btnh = max(ext(:,4)) + 6;
btnsize = repmat([btnw btnh],Nbtns,1);
btnpos = [(0:Nbtns-1).'*(btnw+sep)+sep repmat(sep,Nbtns,1)];

set(hBtn,{'Position'},mat2cell([btnpos btnsize],ones(Nbtns,1),4));
btnpanelsiz = [Nbtns*(btnw+sep)+sep btnh+2*sep];
set(hBtnPanel,'Position',[0 0 btnpanelsiz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output handle struct
handles.fig = fig;
handles.panels = hPanels;
handles.ctrls = hCtrls;
handles.btnpanel = hBtnPanel;
handles.btns = hBtn;

% output positioning info
sinfo.map = map; % Formats grid mapping
sinfo.istext = istext; % true if static text
sinfo.pnminsiz = pnpos(:,[3 4]); % minimum panel size
sinfo.ctrlpos = ctrlpos; % control position w/in panel
sinfo.labelpos = labelpos; % label positions w/in panel
sinfo.labelalign = labelalign; % label alignments
sinfo.btnpnsiz = btnpanelsiz; % button panel size

sinfo.w_max = max(wcell,[],1);
sinfo.w_min = min(wcell,[],1);
sinfo.h_max = max(hcell,[],2);
sinfo.h_min = min(hcell,[],2);

sinfo.w_delta = sinfo.w_max-sinfo.w_min;
sinfo.h_delta = sinfo.h_max-sinfo.h_min;

% total non-adjustable width & height over all controls
sinfo.w_totalnonadj = sum(sinfo.w_max(~autowidth))+sum(sinfo.w_min(autowidth));
sinfo.h_totalnonadj = sum(sinfo.h_max(~autoheight))+sum(sinfo.h_min(autoheight));

% autosize related fields
sinfo.autosize = logical(autosize); % 1 to autosize once, 2 to resize with window
sinfo.autowidth = find(autowidth);
sinfo.autoheight = find(autoheight);

figw = max(sum(sinfo.w_max)+sep*(dim(2)+1),sinfo.btnpnsiz(1)+2*sep);
figh = sum(sinfo.h_max)+sep*(dim(1)+1)+sinfo.btnpnsiz(2);
sinfo.figsizmin = min([figw figh], Options.DialogSizeMinimum);
                
% Set default figure size
Nauto = [numel(sinfo.autowidth) numel(sinfo.autoheight)];
figpos = get(fig,'Position');
if Nauto(1)==0 % all fixed columns
   figpos(3) = sinfo.figsizmin(1);
else
   figpos(3) = max(figpos(3),sinfo.figsizmin(1)+59*Nauto(1));
end
if Nauto(2)==0 % all fixed rows
   figpos(4) = sinfo.figsizmin(2);
else
   figpos(4) = max(figpos(4),sinfo.figsizmin(2)+19*Nauto(2));
end
set(fig,'Position',figpos);
movegui(fig,'onscreen');

end

function resizegui(handles,sinfo,Options)
% This function places all controls in proper place.
% Must be called before the GUI is made visible as buildgui function
% just creates uicontrols and do not place them in proper places.

% KNOWN ISSUE: Resize event does not fire all the time

sep = Options.Sep; % separation between controls

% get current figure size
figPos = get(handles.fig,'Position');
figSize = figPos(3:4);

% force figure size to be larger than the minimum size allowed
idx = figSize<sinfo.figsizmin;
if any(idx)
   figSize(idx) = sinfo.figsizmin(idx);
   figPos([3 4]) = figSize;
   set(handles.fig,'Position',figPos);
   movegui(handles.fig,'onscreen');
end

% Place the button panel (lower right hand corner)
pos = [figSize(1)-sinfo.btnpnsiz(1)-sep sep sinfo.btnpnsiz];
set(handles.btnpanel,'Position',pos);

dim = size(sinfo.map);
num = size(handles.panels,1);

% Determine the column widths
w_total = figSize(1)-sep*(dim(2)+1); % sans spacers
autowidth = sinfo.autowidth;
totalnonadj = sinfo.w_totalnonadj;
w_adj = (w_total-totalnonadj)/numel(autowidth);
widths = sinfo.w_max;
widths(autowidth) = sinfo.w_min(autowidth) + w_adj;

% make sure all the fixed width components fit in the auto-width columns
idx = widths(autowidth)<sinfo.w_max(autowidth);
while any(idx)
   % fix the row with maximum difference b/w max and min
   idx = autowidth(idx); % idx now contains actual row index
   [~,I] = max(sinfo.w_delta(idx)); % max diff
   idx = idx(I); % narrow down to 1
   autowidth = setdiff(autowidth,idx); % remove it from autowidth list
   widths(idx) = sinfo.w_max(idx); % fix the width to its max
   totalnonadj = totalnonadj + sinfo.w_max(idx) - sinfo.w_min(idx);
   
   % recompute the adjustable heights
   w_adj = (w_total-totalnonadj)/numel(autowidth);
   widths(autowidth) = sinfo.w_min(autowidth) + w_adj;
   idx = widths(autowidth)<sinfo.w_max(autowidth);
end

% Determine the grid's x-coordinates
grid_x0 = cumsum([0 widths(1:end-1)]+sep);
grid_x1 = grid_x0 + widths;

% Set heights of autoadjusted text controls
h_min = sinfo.h_min;
heights = sinfo.h_max;
idx = sinfo.istext & sinfo.autosize(:,1);
for m = find(idx).'
   % get the grid cell position
   [i,j] = find(sinfo.map==m);
   x = min(grid_x0(j));
   pwidth = max(grid_x1(j))-x;
   ppos = [x 0 pwidth 1]; % panel position, only set width
   
   set(handles.panels(m,1),'Position',ppos);
   
   % use textwrap to obtain an initial wrapped candidate
   str = get(handles.ctrls(m,1),'String');
   msg = textwrap(handles.ctrls(m,1),{str});
   
   % often too wide, so reduce words on every line.
   % Known Bug: if there is an explicit newline in the string, there
   % is a possible for it to be ignored as a result of the code below
   Nmsgs = numel(msg);
   str = '';
   for k = 1:Nmsgs
      if k==1
         str = msg{1};
      else
         str = sprintf('%s\n%s',str,msg{k});
      end
      set(handles.ctrls(m,2),'String',str);
      ext = get(handles.ctrls(m,2),'Extent');
      while ext(3)>1
         % send last word off to the next line
         [idx,word] = regexp(str,'\s(\S+)$','start','tokens','once');
         str(idx:end) = [];
         set(handles.ctrls(m,2),'String',str);
         ext = get(handles.ctrls(m,2),'Extent');
         if k==Nmsgs && numel(msg)==Nmsgs
            msg(k+1) = word;
         else
            msg{k+1} = sprintf('%s %s',word{1},msg{k+1});
         end
      end
   end
   
   % set height
   i = unique(i);
   heights(i) = max(heights(i),ext(4)/numel(i));
   h_min(i) = min(h_min(i),ext(4)/numel(i));
end

% Determine the automatically adjusted row heights
h_total = figSize(2)-sep*(dim(1)+1)-sinfo.btnpnsiz(2);
autoheight = sinfo.autoheight;
totalnonadj = sinfo.h_totalnonadj;
h_adj = (h_total-totalnonadj)/numel(autoheight);
heights(autoheight) = h_min(autoheight) + h_adj;

idx = heights(autoheight)<sinfo.h_max(autoheight);
while any(idx)
   % fix the row with maximum difference b/w max and min
   idx = autoheight(idx); % idx now contains actual row index
   [~,I] = max(sinfo.h_delta(idx));
   idx = idx(I);
   autoheight = setdiff(autoheight,idx);
   heights(idx) = sinfo.h_max(idx);
   totalnonadj = totalnonadj + sinfo.h_max(idx) - sinfo.h_min(idx);
   
   % recompute adjustable heights
   h_adj = (h_total-totalnonadj)/numel(autoheight);
   heights(autoheight) = sinfo.h_min(autoheight) + h_adj;
   idx = heights(autoheight)<sinfo.h_max(autoheight);
end

% grid y-coordinates
grid_y0 = flipud(cumsum([0;flipud(heights(2:end))]+sep)) + sinfo.btnpnsiz(2);
yoffset = figSize(2) - sep - grid_y0(1) - heights(1); % aligned to top edge
grid_y0(:) = grid_y0 + yoffset;
grid_y1 = grid_y0 + heights;

% obtain position of each control panel
for m = 1:num
   % get the grid cell position
   [i,j] = find(sinfo.map==m);
   x = min(grid_x0(j));
   y = min(grid_y0(i));
   ytop = max(grid_y1(i));
   pwidth = max(grid_x1(j))-x;
   ppos = [x y pwidth ytop-y]; % panel position

   % if fixed height, align the panel to the upper edge of the grid cell
   if ~sinfo.autosize(m,2)
      ppos(2) = ytop-sinfo.pnminsiz(m,2);
      ppos(4) = sinfo.pnminsiz(m,2);
   end
   
   set(handles.panels(m,1),'Position',ppos);
   
   if ~sinfo.istext(m) && any(sinfo.autosize(m,:)) % for all other control types

      cpos = {sinfo.ctrlpos(m,:);sinfo.labelpos(m,:,1);sinfo.labelpos(m,:,2)};
      
      % adjust the width of the control
      if sinfo.autosize(m,1) % auto width
         
         % determine the control width
         dw = ppos(3)-sinfo.pnminsiz(m,1);
         wd = cpos{1}(3) + dw; % new control width
         cpos{1}(3) = wd;
   
         x0 = sinfo.ctrlpos(m,1);
         if ~sinfo.labelalign(m,1) % prompt label at top
            switch sinfo.labelalign(m,2)
               case 2 % right
                  cpos{2}(1) = x0 + wd;
               case 1 % center
                  cpos{2}(1) = x0 + wd/2;
            end
         end
         
         if sinfo.labelalign(m,3) % unit label to right
            cpos{3}(1) = x0 + wd;
         else % unit label at bottom
            switch sinfo.labelalign(m,4)
               case 2 % right
                  cpos{3}(1) = x0 + wd;
               case 1 % center
                  cpos{3}(1) = x0 + wd/2;
            end
         end
      end
      
      if sinfo.autosize(m,2) % auto height
         dh = ppos(4)-sinfo.pnminsiz(m,2);
         ht = cpos{1}(4) + dh;
         cpos{1}(4) = ht;
         
         y0 = sinfo.ctrlpos(m,2);
         if sinfo.labelalign(m,1) % prompt label to left
            switch sinfo.labelalign(m,2)
               case 2 % top
                  cpos{2}(2) = y0 + ht;
               case 1 % middle
                  cpos{2}(2) = y0 + ht/2;
            end
         else % prompt label at top
            cpos{2}(2) = y0 + ht;
         end
         
         if sinfo.labelalign(m,3) % units label to right
            switch sinfo.labelalign(m,4)
               case 2 % top
                  cpos{3}(2) = ht;
               case 1 % middle
                  cpos{3}(2) = ht/2;
            end
         end
      end
      
      set(handles.ctrls(m,:),{'Position'},cpos);
   end
end
end


% Copyright (c) 2009-2015, Takeshi Ikuma
% Copyright (c) 2010, Luke Reisner
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%   * Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%   * Neither the names of its contributors may be used to endorse or
%     promote products derived from this software without specific prior
%     written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
