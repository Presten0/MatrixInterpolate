using System;
using UnityEditor;
using UnityEngine;
using Vectors;

[ExecuteAlways]
[RequireComponent(typeof(VectorRenderer))]
public class ModelScript : MonoBehaviour {
    
    [NonSerialized] 
    private VectorRenderer vectors;

    [SerializeField, Range(0f, 1f)] 
    private float Time;
    
    [SerializeField]
    public Matrix4x4 matrixA = Matrix4x4.identity;
    [SerializeField]
    public Matrix4x4 matrixB = Matrix4x4.identity;
    
    private Matrix4x4 matrixT = Matrix4x4.identity;

    void OnEnable() {
        vectors = GetComponent<VectorRenderer>();
    }

    void Update()
    {
        Vector3 temp1 = matrixA.GetColumn(0);
        Vector3 temp2 = matrixA.GetColumn(1);
        Vector3 temp3 = matrixA.GetColumn(2);
        
        Vector3 temp4 = matrixB.GetColumn(0);
        Vector3 temp5 = matrixB.GetColumn(1);
        Vector3 temp6 = matrixB.GetColumn(2);
        
        Vector3 temp7 = matrixT.GetColumn(0);
        Vector3 temp8 = matrixT.GetColumn(1);
        Vector3 temp9 = matrixT.GetColumn(2);
        
        Debug.Log(GetScaleVector(matrixT));

        var interpolatedTranslation = (1.0f - Time) * GetTranslationVector(matrixA) + Time * GetTranslationVector(matrixB);
        SetTranslationVector(ref matrixT, interpolatedTranslation);
        var interpolatedScale = (1.0f - Time) * GetScaleVector(matrixA) + Time * GetScaleVector(matrixB);
        //SetScale(ref matrixT, interpolatedScale);
        
        using (vectors.Begin()) {
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+temp1, Color.red);
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+temp2, Color.green);
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+temp3, Color.blue);
            
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+temp4, Color.red);
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+temp5, Color.green);
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+temp6, Color.blue);
            
            vectors.Draw(GetTranslationVector(matrixT), GetTranslationVector(matrixT)+temp7, Color.red);
            vectors.Draw(GetTranslationVector(matrixT), GetTranslationVector(matrixT)+temp8, Color.green);
            vectors.Draw(GetTranslationVector(matrixT), GetTranslationVector(matrixT)+temp9, Color.blue);
        }
    }

    public Vector3 GetTranslationVector(Matrix4x4 matrix)
    {
        return matrix.GetColumn(3);
    }

    public void SetTranslationVector(ref Matrix4x4 matrix, Vector3 translationVector)
    {
        matrix.SetColumn(3, translationVector);
    }

    public Vector3 GetScaleVector(Matrix4x4 matrix)
    {
        var scaleX = MathF.Sqrt(matrix.m00 * matrix.m00 + matrix.m01 * matrix.m01 + matrix.m02 * matrix.m02);
        var scaleY = MathF.Sqrt(matrix.m10 * matrix.m10 + matrix.m11 * matrix.m11 + matrix.m12 * matrix.m12);
        var scaleZ = MathF.Sqrt(matrix.m20 * matrix.m20 + matrix.m21 * matrix.m21 + matrix.m22 * matrix.m22);

        return new Vector3(scaleX, scaleY, scaleZ);
    }

    public void SetScale(ref Matrix4x4 matrix, Vector3 scaleVector)
    {
        var column1 = matrix.GetColumn(0);
        var column2 = matrix.GetColumn(1);
        var column3 = matrix.GetColumn(2);

        var scalexFactor = column1.x * scaleVector.x;
        var scaleyFactor = column2.y * scaleVector.y;
        var scalezFactor = column3.z * scaleVector.z;

        matrix.m00 = scalexFactor;
        matrix.m01 = scaleyFactor;
        matrix.m02 = scalezFactor;

    }
}

[CustomEditor(typeof(ModelScript))]
public class ExampleGUI : Editor {
    void OnSceneGUI() {
        var ex = target as ModelScript;
        if (ex == null) return;

        EditorGUI.BeginChangeCheck();
        var a = Handles.PositionHandle(ex.GetTranslationVector(ex.matrixA), Quaternion.identity); //Handle at translation-vector
        var b = Handles.PositionHandle(ex.GetTranslationVector(ex.matrixB), Quaternion.identity); //Handle at translation-vector

        if (EditorGUI.EndChangeCheck()) {
            Undo.RecordObject(target, "Vector Positions");
            ex.SetTranslationVector(ref ex.matrixA, a); //Set translation-vector to handles position
            ex.SetTranslationVector(ref ex.matrixB, b); //Set translation-vector to handles position
            EditorUtility.SetDirty(target);
        }
    }
}
