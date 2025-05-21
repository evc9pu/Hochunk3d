Function Read_Array,FileName,Max_Wave,Min_Wave,Max_Flux,Min_Flux

Fluxes = Replicate({model_output, Wave:FltArr(200), Flux:DblArr(200)},10,3)
OpenR,1,FileName

Max_Wave = 0
Min_Wave = 1000
Max_Flux = 0
Min_Flux = 1

Header = ''
ReadF,1,Header
For Aperture = 0,2 Do Begin
  For Angle = 0,9 Do Begin
    For Result = 0,199 Do Begin
      ReadF,1,Wave,Flux
      Fluxes[Angle,Aperture].Wave[Result] = Wave
      Fluxes[Angle,Aperture].Flux[Result] = Flux
      If (Flux GT Max_Flux) Then Max_Flux = Flux
      If ((Flux NE 0) And (Flux LT Min_Flux)) Then Min_Flux = Flux
    EndFor
  EndFor
EndFor
Max_Wave = Max(Fluxes.Wave[199])
Min_Wave = Min(Fluxes.Wave[0])

Close,1
Return,Fluxes

End
