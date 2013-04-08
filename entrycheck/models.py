from django.db import models
from django import forms
#from widgets import TextArea 
from django.forms import ValidationError
#import json

# Create your models here.
class EntryForm(forms.Form):
    import helptexts
    out=""
    stuff={}
    t1=[""]
    t2=[""]
    t3=[""]
    name = forms.CharField(max_length=200, initial='', help_text=helptexts.helptxt['name'])
    env=forms.BooleanField(required=False)
    addHs=forms.BooleanField(required=False)
    reaction =  forms.CharField(max_length=400, widget=forms.Textarea( attrs={'rows':'2', 'cols': '10'}), initial='CCO', help_text=helptexts.helptxt['reaction'])
    test1 = forms.CharField(max_length=400, widget=forms.Textarea( attrs={'rows':'2', 'cols': '10'}),  initial='', help_text=helptexts.helptxt['test'])
    test2 = forms.CharField(max_length=400, widget=forms.Textarea( attrs={'rows':'2', 'cols': '10'}), initial='', help_text=helptexts.helptxt['test'])
    test3 = forms.CharField(max_length=400, widget=forms.Textarea( attrs={'rows':'2', 'cols': '10'}), initial='', help_text=helptexts.helptxt['test'])
    products = forms.CharField(max_length=200, widget=forms.Textarea( attrs={'rows':'2', 'cols': '10'}), initial='', help_text=helptexts.helptxt['products'])
    conditions =  forms.CharField(max_length=400,widget=forms.Textarea( attrs={'rows':'2', 'cols': '10'}), initial='', help_text=helptexts.helptxt['conditions'])
    synthons = forms.IntegerField(help_text=helptexts.helptxt['synthons'])
    doi =  forms.CharField(max_length=400, initial='', help_text=helptexts.helptxt['doi'])
    done=forms.BooleanField(required=True)
        
    def clean(self):
       cleaned_data = self.cleaned_data
       if not self.seriouslyValidate(cleaned_data):
          raise forms.ValidationError("invalid Entry")
       return cleaned_data # Never forget this otherwise you will end up with empty self.cleaned_data.get( in the views.py.)

    def seriouslyValidate(self,cleaned_data):
        import chemtests as check
        result = True
        if not check.checkName(cleaned_data.get('name')):
            raise forms.ValidationError("invalid Name")
            result = False
            
        reaction1ok, scriptout, stuff1=check.checkReaction(cleaned_data.get('reaction'), cleaned_data.get('test1'), cleaned_data.get('addHs'), cleaned_data.get('env'))
        self.out= ' <b>test1:</b> ' + scriptout
        self.stuff['test1']=stuff1
        if stuff1 is not None and stuff1 != {}:
            #t1 etc should become lists, which should be printed since the same reaction can occur on multiple spots
            #self.t1=[]
            self.t1=[]
            for i in range(len(stuff1['products'])): 
                self.t1.append('.'.join(stuff1['reactants']) + '>>' + stuff1['products'][i])#+ pt1[0]#'.'.join(stuff1['products'])
            print self.t1
        reaction2ok, scriptout2, stuff2=check.checkReaction(cleaned_data.get('reaction'), cleaned_data.get('test2'), cleaned_data.get('addHs'), cleaned_data.get('env'))
        self.out=self.out +' <b>test2:</b> ' + scriptout2
        self.stuff['test2']=stuff2
        if stuff2 is not None and stuff2 != {}:
            self.t2=[]
            for i in range(len(stuff2['products'])): 
                self.t2.append('.'.join(stuff2['reactants']) + '>>' + stuff2['products'][i])
            print self.t2
        reaction3ok, scriptout3, stuff3=check.checkReaction(cleaned_data.get('reaction'), cleaned_data.get('test3'), cleaned_data.get('addHs'), cleaned_data.get('env'))
        self.out=self.out +' <b>test3:</b> ' + scriptout3
        self.stuff['test3']=stuff3
        if stuff3 is not None and stuff3 != {}:
            self.t3=[]
            for i in range(len(stuff3['products'])): 
                self.t3.append('.'.join(stuff3['reactants']) + '>>' + stuff3['products'][i]) 
            print self.t3
        if not reaction1ok:     
            #print scriptout
            result = False
            raise forms.ValidationError("invalid Reaction with test 1")
        
        if not reaction2ok:     
            #print scriptout
            result = False
            raise forms.ValidationError(" invalid Reaction with test 2 ")
        
        if not reaction3ok:     
            #print scriptout
            result = False
            raise forms.ValidationError(" invalid Reaction with test 3 ")
               
               
               
        if not check.checkProducts(cleaned_data.get('products')):
            raise forms.ValidationError("invalid Products")
            result = False
        if not check.checkConditions(cleaned_data.get('conditions')):
            raise forms.ValidationError("invalid conditions")
            result = False
        if not check.checkSynthons(cleaned_data.get('synthons')):
            raise forms.ValidationError("invalid synthons")
            result = False
        if not check.checkDoi(cleaned_data.get('doi')):
            raise forms.ValidationError("invalid doi")
            result = False
        
        return result
    
    def getOut(self):
        return self.out
        
    def getStuff(self):
        return self.stuff
    def getFullReactions(self):
           return [self.t1,self.t2,self.t3]
        