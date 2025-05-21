Pro Handle_Models_Event,Event

Common Models_Block,Models_1,Models_2,Models_3,Models_4,Models_5,Active_Model,Which_Aper,Which_Angle
Common Plot_Param,Max_Wave,Min_Wave,Max_Flux,Min_Flux,Which_Model
Common User_Wid,button_aperture,button_angle

Widget_Control,Event.ID,Get_UValue=uvalue
Case UValue Of
  'quit' : Widget_Control,Event.TOP,/Destroy
  'plot' : Begin
             Plot_OO,[Min_Wave,Max_Wave],[Min_Flux,Max_Flux],/NoData
             For I=0,4 Do Begin
               For J=0,9 Do Begin
                 For K=0,4 Do Begin
                   If (Which_Model[K]) Then Begin
                     Case K Of
                       0 : Models = Models_1
                       1 : Models = Models_2
                       2 : Models = Models_3
                       3 : Models = Models_4
                       4 : Models = Models_5
                     EndCase
                     If (Which_Aper[K,I] And Which_Angle[K,J]) Then Begin
                       Index = Where(Models[J,I].Flux)
                       OPlot,Models[J,I].Wave[Index],Models[J,I].Flux[Index]
                     EndIf
                   EndIf
                 EndFor
               EndFor
             EndFor
           End
  'plot_c' : Begin
               Plot,[-1,1],[-1,1],/NoData
               For I=0,4 Do Begin
                 For J=0,9 Do Begin
                   For K=0,4 Do Begin
                     If (Which_Model[K]) Then Begin
                       Case K Of
                         0 : Models = Models_1
                         1 : Models = Models_2
                         2 : Models = Models_3
                         3 : Models = Models_4
                         4 : Models = Models_5
                       EndCase
                       If (Which_Aper[K,I] And Which_Angle[K,J]) Then Begin
                         Dummy = Min(Abs(Models[J,I].Wave-60.),L_60)
                         Dummy = Min(Abs(Models[J,I].Wave-100.),L_100)
                         Dummy = Min(Abs(Models[J,I].Wave-170.),L_170)
                         C1 = ALog10(Models[J,I].Flux[L_100]/Models[J,I].Flux[L_60])
                         C2 = ALog10(Models[J,I].Flux[L_170]/Models[J,I].Flux[L_100])
                         PlotS,C2,C1,PSym=3
                       EndIf
                     EndIf
                   EndFor
                 EndFor
               EndFor
             End
  'Sel_Angle' : Which_Angle[Active_Model,Event.Value] = Event.Select
  'Sel_Aper' : Which_Aper[Active_Model,Event.Value] = Event.Select
  'Sel_Plot' : Which_Model[Event.Value] = Event.Select
EndCase

Return
End

Pro Text_Sel,Wid_Str

Common Models_Block,Models_1,Models_2,Models_3,Models_4,Models_5,Active_Model,Which_Aper,Which_Angle
Common Plot_Param,Max_Wave,Min_Wave,Max_Flux,Min_Flux,Which_Model

If (Tag_Names(Wid_Str,/STRUCTURE_NAME) NE 'WIDGET_CONTEXT') Then Return
Result = Dialog_PickFile(Dialog_Parent=base,Filter='*.dat',/MUST_EXIST)
If (Result Eq '') Then Return
Widget_Control,Wid_Str.ID,Set_Value=Result
Widget_Control,Wid_Str.ID,Get_UValue=uvalue
Models = Read_Array(Result,Max_Wave,Min_Wave,Max_Flux,Min_Flux)
Min_Flux=Max_flux*1.e-5
Case UValue Of
  'mod_1' : Models_1 = Models
  'mod_2' : Models_2 = Models
  'mod_3' : Models_3 = Models
  'mod_4' : Models_4 = Models
  'mod_5' : Models_5 = Models
EndCase

Return
End

Function Sel_Mod,Which_Button

Common Models_Block,Models_1,Models_2,Models_3,Models_4,Models_5,Active_Model,Which_Aper,Which_Angle
Common User_Wid,button_aperture,button_angle

Active_Model = Which_Button.Value
Widget_Control,button_aperture,Set_Value=Which_Aper[Active_Model,*]
Widget_Control,button_angle,Set_Value=Which_Angle[Active_Model,*]

Return,1
End

Pro Handle_Models

Common Models_Block,Models_1,Models_2,Models_3,Models_4,Models_5,Active_Model,Which_Aper,Which_Angle
Common User_Wid,button_aperture,button_angle
Common Plot_Param,Max_Wave,Min_Wave,Max_Flux,Min_Flux,Which_Model

Device,Decomposed=0
Which_Aper = Replicate(0B,5,5)
Which_Angle = Replicate(0B,5,10)
Which_Model = Replicate(0B,5)
Angles = ['87','81','76','70','63','57','49','41','32','18']
base = Widget_Base(Row=3)
plot_models = Widget_Draw(base,XSize=640,Ysize=512)
base2 = Widget_Base(base,Row=2)
base3 = Widget_Base(base,/Column)
base4 = Widget_Base(base,Column=2)
button_aperture = CW_BGroup(base2,['1','2','3','4','5'],/Column,/NonExclusive,Label_top='Aperture',$
                  UValue='Sel_Aper',/Frame)
button_angle = CW_BGroup(base2,Angles,/Column,/NonExclusive,Label_top='Angle',/Frame,$
               UValue='Sel_Angle')
button_text_1 = Widget_Text(base3, VALUE='File 1 ...',/Context_Events,Event_Pro='Text_Sel',$
                /All_Events,XSize=100,UValue='mod_1')
button_text_2 = Widget_Text(base3, VALUE='File 2 ...',/Context_Events,Event_Pro='Text_Sel',$
                /All_Events,XSize=100,UValue='mod_2')
button_text_3 = Widget_Text(base3, VALUE='File 3 ...',/Context_Events,Event_Pro='Text_Sel',$
                /All_Events,XSize=100,UValue='mod_3')
button_text_4 = Widget_Text(base3, VALUE='File 4 ...',/Context_Events,Event_Pro='Text_Sel',$
                /All_Events,XSize=100,UValue='mod_4')
button_text_5 = Widget_Text(base3, VALUE='File 5 ...',/Context_Events,Event_Pro='Text_Sel',$
                /All_Events,XSize=100,UValue='mod_5')
button_sm = CW_BGroup(base4,Replicate('',5),/Column,/Exclusive,Label_top='Sel',$
                  Event_Func='Sel_Mod',Set_Value=0,/Frame,YSize=150)
Active_Model = 0
button_pm = CW_BGroup(base4,Replicate('',5),/Column,/NonExclusive,Label_top='Plot',$
                  UValue='Sel_Plot',/Frame)
button_quit = Widget_Button(base2, VALUE='End application',UValue='quit')
button_plot = Widget_Button(base, VALUE='Plot selected models',UValue='plot')
;button_plot_c = Widget_Button(base, VALUE='Plot colors of selected models',UValue='plot_c')

Widget_Control, base, /Realize
XManager,'Handle_Models',base

Return
End
