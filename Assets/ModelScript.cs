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
    
    [SerializeField,HideInInspector]internal Matrix4x4 matrixA = Matrix4x4.identity;
    [SerializeField,HideInInspector]internal Matrix4x4 matrixB = Matrix4x4.identity;
    [SerializeField,HideInInspector]internal Matrix4x4 matrixI = Matrix4x4.identity;

    void OnEnable() {
        vectors = GetComponent<VectorRenderer>();
    }

    void Update()
    {
        //matrixA = Matrix4x4.identity;
        //matrixB = Matrix4x4.identity;
        
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
        matrix.m03 = translationVector.x;
        matrix.m13 = translationVector.y;
        matrix.m23 = translationVector.z;
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

    public Quaternion GetRotation(Matrix4x4 matrix)
    {
        var cosX = Vector3.Dot(GetNormalizedVector(matrix.GetColumn(0)),new Vector3(1, 0, 0));
        var cosY = Vector3.Dot(GetNormalizedVector(matrix.GetColumn(1)),new Vector3(0, 1, 0));
        var cosZ = Vector3.Dot(GetNormalizedVector(matrix.GetColumn(2)),new Vector3(0, 0, 1));

        var smallCross = Vector3.zero;
        var sinHalf = 0.0f;
        var cosHalf = 0.0f;
        
        if (cosX < cosY && cosX < cosZ)
        {
            smallCross = Vector3.Cross(GetNormalizedVector(matrix.GetColumn(0)),new Vector3(1, 0, 0));
            sinHalf= MathF.Sqrt((1 - cosX) / 2);
            cosHalf= MathF.Sqrt((1 + cosX) / 2);

        } else if (cosY < cosX && cosY < cosZ)
        {
            smallCross = Vector3.Cross(GetNormalizedVector(matrix.GetColumn(1)),new Vector3(0, 1, 0));
            sinHalf= MathF.Sqrt((1 - cosY) / 2);
            cosHalf= MathF.Sqrt((1 + cosY) / 2);
            
        } else if (cosZ < cosX && cosZ < cosY)
        {
            smallCross = Vector3.Cross(GetNormalizedVector(matrix.GetColumn(2)),new Vector3(0, 0, 1));
            sinHalf= MathF.Sqrt((1 - cosZ) / 2);
            cosHalf= MathF.Sqrt((1 + cosZ) / 2);
            
        } else
        {
            return Quaternion.identity;
        }

        var quat = smallCross * sinHalf;

        return new Quaternion(quat.x, quat.y, quat.z, cosHalf);
    }

    public Quaternion GetInterpolatedQuaternion(Quaternion quatA, Quaternion quatB)
    {
        var quatC = new Quaternion(quatA.x * quatB.x, quatA.y * quatB.y, quatA.z * quatB.z, -quatA.w * quatB.w);
        
        var lerpedAngle = MathF.Acos((1.0f - Time) * quatA.w + Time * quatB.w);
    }

    public Vector3 GetNormalizedVector(Vector3 vector)
    {
        var magnitude = MathF.Sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
        return vector / magnitude;
    }
}

[CustomEditor(typeof(ModelScript))]
public class ModelGUI : Editor {
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
    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();
        
        var ex = target as ModelScript;
        if (ex == null) return;

        EditorGUI.BeginChangeCheck();

        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Matrix A");
        EditorGUILayout.BeginVertical();

        var resultA = Matrix4x4.identity;
        for (var i = 0; i < 4; i++)
        {
            EditorGUILayout.BeginHorizontal();
            for (var j = 0; j < 4; j++)
            {
                resultA[i, j] = EditorGUILayout.FloatField((ex.matrixA[i, j]));
            }
            EditorGUILayout.EndHorizontal();
        }
        
        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();
        
        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Matrix B");
        EditorGUILayout.BeginVertical();

        var resultB = Matrix4x4.identity;
        for (var i = 0; i < 4; i++)
        {
            EditorGUILayout.BeginHorizontal();
            for (var j = 0; j < 4; j++)
            {
                resultB[i, j] = EditorGUILayout.FloatField((ex.matrixB[i, j]));
            }
            EditorGUILayout.EndHorizontal();
        }
        
        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();
        
        if (EditorGUI.EndChangeCheck())
        {
            Undo.RecordObject(ex, "Changed matrix");
            ex.matrixA = resultA;
            ex.matrixB = resultB;
            EditorUtility.SetDirty(ex);
        }
    }
}


