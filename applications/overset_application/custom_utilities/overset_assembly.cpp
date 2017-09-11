

class OversetAssembly
{
public:
    OversetAssembly(ModelPart &, Communicator &);

    ~OversetAssembly();

    //volume to face
    //face to volume
    GenerateAdditionalConnectivity()
    {
        

    }

    //crd, cnn
    GenerateCrdCnn();

    //searcher
    GenerateSearcher();

    SearchHinges();

    //
    GetHingesValues();

private:
    ModelPart & mrModelPart;


    
}