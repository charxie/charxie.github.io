//periodic table creator

 function moveover(txt){ 
    window.status = txt; 
 } 
 function fillitin(Name, AtomicNumber,AtomicWeight, Shells, FillingOrbital, MeltingPoint, FreezingPoint){ 
    moveover(Name); 
    document.PeriodicTable.Name.value=Name; 
    document.PeriodicTable.AtomicNumber.value=AtomicNumber; 
    document.PeriodicTable.AtomicWeight.value=AtomicWeight; 
    document.PeriodicTable.Shells.value=Shells; 
    document.PeriodicTable.FillingOrbital.value=FillingOrbital; 
    document.PeriodicTable.MeltingPoint.value=MeltingPoint; 
    document.PeriodicTable.FreezingPoint.value=FreezingPoint; 
 } 
