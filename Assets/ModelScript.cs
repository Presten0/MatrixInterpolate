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
    
    private Matrix4x4 matrixI = Matrix4x4.identity;

    void OnEnable() {
        vectors = GetComponent<VectorRenderer>();
    }

    void Update()
    {
        var interpolatedTranslation = (1.0f - Time) * GetTranslationVector(matrixA) + Time * GetTranslationVector(matrixB);
        SetTranslationVector(ref matrixI, interpolatedTranslation);
        /*var interpolatedScale = (1.0f - Time) * GetScaleVector(matrixA) + Time * GetScaleVector(matrixB);
        SetScale(ref matrixT, interpolatedScale);*/
        
        Vector3 AX = matrixA.GetColumn(0);
        Vector3 AY = matrixA.GetColumn(1);
        Vector3 AZ = matrixA.GetColumn(2);
        
        Vector3 BX = matrixB.GetColumn(0);
        Vector3 BY = matrixB.GetColumn(1);
        Vector3 BZ = matrixB.GetColumn(2);
        
        Vector3 IX = matrixI.GetColumn(0);
        Vector3 IY = matrixI.GetColumn(1);
        Vector3 IZ = matrixI.GetColumn(2);
        
        Debug.Log(GetInterpolatedScaleMagnitudes(matrixA, matrixB));

        using (vectors.Begin()) {
            
            //Visualize matrix A
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+AX, Color.red);
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+AY, Color.green);
            vectors.Draw(GetTranslationVector(matrixA), GetTranslationVector(matrixA)+AZ, Color.blue);
            
            //Visualize matrix B
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+BX, Color.red);
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+BY, Color.green);
            vectors.Draw(GetTranslationVector(matrixB), GetTranslationVector(matrixB)+BZ, Color.blue);
            
            //Visualize box from interpolated matrix
            vectors.Draw(GetTranslationVector(matrixI), GetTranslationVector(matrixI)+IX, Color.red);
            vectors.Draw(GetTranslationVector(matrixI)+IY, GetTranslationVector(matrixI)+IX+IY, Color.red);
            vectors.Draw(GetTranslationVector(matrixI)+IZ, GetTranslationVector(matrixI)+IX+IZ, Color.red);
            vectors.Draw(GetTranslationVector(matrixI)+IY+IZ, GetTranslationVector(matrixI)+IX+IY+IZ, Color.red);
            vectors.Draw(GetTranslationVector(matrixI), GetTranslationVector(matrixI)+IY, Color.green);
            vectors.Draw(GetTranslationVector(matrixI)+IX, GetTranslationVector(matrixI)+IY+IX, Color.green);
            vectors.Draw(GetTranslationVector(matrixI)+IZ, GetTranslationVector(matrixI)+IY+IZ, Color.green);
            vectors.Draw(GetTranslationVector(matrixI)+IX+IZ, GetTranslationVector(matrixI)+IY+IX+IZ, Color.green);
            vectors.Draw(GetTranslationVector(matrixI), GetTranslationVector(matrixI)+IZ, Color.blue);
            vectors.Draw(GetTranslationVector(matrixI)+IX, GetTranslationVector(matrixI)+IZ+IX, Color.blue);
            vectors.Draw(GetTranslationVector(matrixI)+IY, GetTranslationVector(matrixI)+IZ+IY, Color.blue);
            vectors.Draw(GetTranslationVector(matrixI)+IX+IY, GetTranslationVector(matrixI)+IZ+IX+IY, Color.blue);
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

    public Vector3 GetInterpolatedScaleMagnitudes(Matrix4x4 matrixA, Matrix4x4 matrixB)
    {
        var scaleAX = MathF.Sqrt(matrixA.m00 * matrixA.m00 + matrixA.m10 * matrixA.m10 + matrixA.m20 * matrixA.m20);
        var scaleAY = MathF.Sqrt(matrixA.m01 * matrixA.m01 + matrixA.m11 * matrixA.m11 + matrixA.m21 * matrixA.m21);
        var scaleAZ = MathF.Sqrt(matrixA.m02 * matrixA.m02 + matrixA.m12 * matrixA.m12 + matrixA.m22 * matrixA.m22);
        
        var scaleBX = MathF.Sqrt(matrixB.m00 * matrixB.m00 + matrixB.m01 * matrixB.m01 + matrixB.m02 * matrixB.m02);
        var scaleBY = MathF.Sqrt(matrixB.m01 * matrixB.m01 + matrixB.m11 * matrixB.m11 + matrixB.m21 * matrixB.m21);
        var scaleBZ = MathF.Sqrt(matrixB.m02 * matrixB.m02 + matrixB.m12 * matrixB.m12 + matrixB.m22 * matrixB.m22);

        var scaleIX = (1.0f - Time) * scaleAX + Time * scaleBX;
        var scaleIY = (1.0f - Time) * scaleAY + Time * scaleBY;
        var scaleIZ = (1.0f - Time) * scaleAZ + Time * scaleBZ;
        
        return new Vector3(scaleIX, scaleIY, scaleIZ);
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
