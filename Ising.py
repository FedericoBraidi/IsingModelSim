import math
import random
import numpy as np
class Ising:
    
    def __init__ (self,temp, N=100):
        
        self.N = N
        self.kb = 1.380649e-23
        self.temp = temp
        self.J = 1
        self.table = 2*(np.random.randint(2,size=(N,N))) - 1
        self.Et = 0
        self.Mt = 0
        self.Ct = 0
        
    def evoluzione (self):
        """
        Fa un passo di evoluzione (1 flip)
        """
        for i in range (self.N**2):
            riga=random.randint(0,self.N-1)
            colonna=random.randint(0,self.N-1)
            delta_en=self.variazione_energia(riga,colonna)
            if(self.flip_accettato(delta_en)):
                self.table[riga,colonna]*=-1
                self.Ct += delta_en**2 + (2*delta_en*self.Et)
                self.Et += delta_en
                self.Mt -= 2*self.table[riga,colonna]
        return self.table
            
    def variazione_energia (self,riga,colonna):
        """
        Si passano come parametri la riga e la colonna 
        dello spin cambiato per semplificare il calcolo
        - en_primi_ è l'energia dovuta ai primi vicini 
            prima e dopo il flip dello spin
        """
        primi = self.table[(riga-1)%self.N,colonna]+self.table[(riga+1)%self.N,colonna]+\
                          self.table[riga,(colonna-1)%self.N]+self.table[riga,(colonna+1)%self.N]
        s = self.table[riga,colonna]
        delta_en = 2*s*primi
        return delta_en
    
    def flip_accettato(self,delta_en):
        """
        Si passa la delta_energia per calcolare la 
        probabilità di accettare la mossa
        """
        accettato = False
        exp = math.exp(-(delta_en/self.temp))
        if(delta_en<0):
            accettato = True
        else:
            x = random.random()
            if(exp>x):
                accettato = True
            else:
                accettato = False
        return accettato
        
        
    def energia(self):
        """
        Calcola l'energia totale del sistema in una data configurazione
        """
        en = 0
        for i in range(len(self.table)):
            for j in range(len(self.table)):
                primi_vicini = self.table[(i-1)%self.N,j]+self.table[(i+1)%self.N,j]+\
                          self.table[i,(j-1)%self.N]+self.table[i,(j+1)%self.N]
                en -= self.table[i,j]*primi_vicini
        return en/4.
    
    def magnetizzazione(self):
        """
        Calcola la magnetizzazione del sistema in una data configurazione
        """
        return np.sum(self.table)
    
    def get_energia(self):
        
        return self.Et
    
    def get_magnetizzazione(self):
        
        return self.Mt
    
    def get_table(self):
        
        return self.table

    def capacita(self):
        
        C = self.energia()**2
        return C
    
    def get_capacita(self):
        
        return self.Ct